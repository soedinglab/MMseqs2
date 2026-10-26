/**
 * @file schema.c
 * @brief Schema management
 */

#include "core/allocator.h"
#include <carquet/carquet.h>
#include "reader/reader_internal.h"
#include "thrift/parquet_types.h"
#include "core/arena.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ============================================================================
 * Schema Creation
 * ============================================================================
 */

/* Initial and growth capacity for schema arrays */
#define SCHEMA_INITIAL_CAPACITY 64
#define SCHEMA_GROWTH_FACTOR 2

carquet_schema_t* carquet_schema_create(carquet_error_t* error) {
    carquet_schema_t* schema = carquet_mem_calloc(1, sizeof(carquet_schema_t));
    if (!schema) {
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate schema");
        return NULL;
    }

    if (carquet_arena_init_size(&schema->arena, 4096) != CARQUET_OK) {
        carquet_mem_free(schema);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate schema arena");
        return NULL;
    }

    /* Allocate initial arrays with malloc (supports realloc for growth) */
    schema->capacity = SCHEMA_INITIAL_CAPACITY;
    schema->num_elements = 1;  /* Root element */

    schema->elements = carquet_mem_calloc(schema->capacity, sizeof(parquet_schema_element_t));
    if (!schema->elements) {
        carquet_arena_destroy(&schema->arena);
        carquet_mem_free(schema);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate schema elements");
        return NULL;
    }

    /* Initialize root element */
    schema->elements[0].name = carquet_arena_strdup(&schema->arena, "schema");
    schema->elements[0].num_children = 0;

    /* Allocate parent index tracking */
    schema->parent_indices = carquet_mem_calloc(schema->capacity, sizeof(int32_t));
    if (!schema->parent_indices) {
        carquet_mem_free(schema->elements);
        carquet_arena_destroy(&schema->arena);
        carquet_mem_free(schema);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate parent indices");
        return NULL;
    }
    schema->parent_indices[0] = -1;  /* Root has no parent */

    /* Allocate leaf tracking arrays with malloc */
    schema->leaf_indices = carquet_mem_calloc(schema->capacity, sizeof(int32_t));
    schema->max_def_levels = carquet_mem_calloc(schema->capacity, sizeof(int16_t));
    schema->max_rep_levels = carquet_mem_calloc(schema->capacity, sizeof(int16_t));
    schema->num_leaves = 0;

    if (!schema->leaf_indices || !schema->max_def_levels || !schema->max_rep_levels) {
        carquet_mem_free(schema->elements);
        carquet_mem_free(schema->parent_indices);
        carquet_mem_free(schema->leaf_indices);
        carquet_mem_free(schema->max_def_levels);
        carquet_mem_free(schema->max_rep_levels);
        carquet_arena_destroy(&schema->arena);
        carquet_mem_free(schema);
        CARQUET_SET_ERROR(error, CARQUET_ERROR_OUT_OF_MEMORY, "Failed to allocate schema leaf arrays");
        return NULL;
    }

    return schema;
}

void carquet_schema_free(carquet_schema_t* schema) {
    if (schema) {
        carquet_mem_free(schema->elements);
        carquet_mem_free(schema->parent_indices);
        carquet_mem_free(schema->leaf_indices);
        carquet_mem_free(schema->max_def_levels);
        carquet_mem_free(schema->max_rep_levels);
        carquet_arena_destroy(&schema->arena);
        carquet_mem_free(schema);
    }
}

/* Helper to grow schema arrays when capacity is reached */
static carquet_status_t schema_ensure_capacity(carquet_schema_t* schema, int32_t required) {
    if (required <= schema->capacity) {
        return CARQUET_OK;
    }

    int32_t new_capacity = schema->capacity;
    while (new_capacity < required) {
        new_capacity *= SCHEMA_GROWTH_FACTOR;
    }

    parquet_schema_element_t* new_elements = carquet_mem_realloc(
        schema->elements, new_capacity * sizeof(parquet_schema_element_t));
    if (!new_elements) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    /* Zero the new portion */
    memset(new_elements + schema->capacity, 0,
           (new_capacity - schema->capacity) * sizeof(parquet_schema_element_t));
    schema->elements = new_elements;

    int32_t* new_parent_indices = carquet_mem_realloc(
        schema->parent_indices, new_capacity * sizeof(int32_t));
    if (!new_parent_indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    schema->parent_indices = new_parent_indices;

    int32_t* new_leaf_indices = carquet_mem_realloc(
        schema->leaf_indices, new_capacity * sizeof(int32_t));
    if (!new_leaf_indices) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    schema->leaf_indices = new_leaf_indices;

    int16_t* new_max_def = carquet_mem_realloc(
        schema->max_def_levels, new_capacity * sizeof(int16_t));
    if (!new_max_def) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    schema->max_def_levels = new_max_def;

    int16_t* new_max_rep = carquet_mem_realloc(
        schema->max_rep_levels, new_capacity * sizeof(int16_t));
    if (!new_max_rep) {
        return CARQUET_ERROR_OUT_OF_MEMORY;
    }
    schema->max_rep_levels = new_max_rep;

    schema->capacity = new_capacity;
    return CARQUET_OK;
}

/* ============================================================================
 * Schema Building
 * ============================================================================
 */

static int32_t decimal_max_precision_for_fixed_len(int32_t type_length) {
    if (type_length <= 0) {
        return 0;
    }

    long double bits = (long double)type_length * 8.0L - 1.0L;
    return (int32_t)floorl(bits * log10l(2.0L));
}

static carquet_status_t validate_column_logical_type(
    carquet_physical_type_t physical_type,
    const carquet_logical_type_t* logical_type,
    int32_t type_length) {

    if (!logical_type || logical_type->id == CARQUET_LOGICAL_UNKNOWN) {
        return CARQUET_OK;
    }

    switch (logical_type->id) {
        case CARQUET_LOGICAL_STRING:
        case CARQUET_LOGICAL_ENUM:
        case CARQUET_LOGICAL_JSON:
        case CARQUET_LOGICAL_BSON:
        case CARQUET_LOGICAL_GEOMETRY:
        case CARQUET_LOGICAL_GEOGRAPHY:
            return physical_type == CARQUET_PHYSICAL_BYTE_ARRAY
                ? CARQUET_OK
                : CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_DATE:
            return physical_type == CARQUET_PHYSICAL_INT32
                ? CARQUET_OK
                : CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_TIME:
            if (logical_type->params.time.unit == CARQUET_TIME_UNIT_MILLIS) {
                return physical_type == CARQUET_PHYSICAL_INT32
                    ? CARQUET_OK
                    : CARQUET_ERROR_INVALID_ARGUMENT;
            }
            if (logical_type->params.time.unit == CARQUET_TIME_UNIT_MICROS ||
                logical_type->params.time.unit == CARQUET_TIME_UNIT_NANOS) {
                return physical_type == CARQUET_PHYSICAL_INT64
                    ? CARQUET_OK
                    : CARQUET_ERROR_INVALID_ARGUMENT;
            }
            return CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_TIMESTAMP:
            return physical_type == CARQUET_PHYSICAL_INT64
                ? CARQUET_OK
                : CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_INTEGER:
            switch (logical_type->params.integer.bit_width) {
                case 8:
                case 16:
                case 32:
                    return physical_type == CARQUET_PHYSICAL_INT32
                        ? CARQUET_OK
                        : CARQUET_ERROR_INVALID_ARGUMENT;
                case 64:
                    return physical_type == CARQUET_PHYSICAL_INT64
                        ? CARQUET_OK
                        : CARQUET_ERROR_INVALID_ARGUMENT;
                default:
                    return CARQUET_ERROR_INVALID_ARGUMENT;
            }

        case CARQUET_LOGICAL_DECIMAL: {
            int32_t precision = logical_type->params.decimal.precision;
            int32_t scale = logical_type->params.decimal.scale;
            if (precision <= 0 || scale < 0 || scale > precision) {
                return CARQUET_ERROR_INVALID_ARGUMENT;
            }

            switch (physical_type) {
                case CARQUET_PHYSICAL_INT32:
                    return precision <= 9 ? CARQUET_OK : CARQUET_ERROR_INVALID_ARGUMENT;
                case CARQUET_PHYSICAL_INT64:
                    return precision <= 18 ? CARQUET_OK : CARQUET_ERROR_INVALID_ARGUMENT;
                case CARQUET_PHYSICAL_BYTE_ARRAY:
                    return CARQUET_OK;
                case CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY:
                    return precision <= decimal_max_precision_for_fixed_len(type_length)
                        ? CARQUET_OK
                        : CARQUET_ERROR_INVALID_ARGUMENT;
                default:
                    return CARQUET_ERROR_INVALID_ARGUMENT;
            }
        }

        case CARQUET_LOGICAL_UUID:
            return physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY &&
                   type_length == 16
                ? CARQUET_OK
                : CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_FLOAT16:
            return physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY &&
                   type_length == 2
                ? CARQUET_OK
                : CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_INTERVAL:
            /* INTERVAL is a 12-byte FIXED_LEN_BYTE_ARRAY (months/days/millis). */
            return physical_type == CARQUET_PHYSICAL_FIXED_LEN_BYTE_ARRAY &&
                   type_length == 12
                ? CARQUET_OK
                : CARQUET_ERROR_INVALID_ARGUMENT;

        case CARQUET_LOGICAL_NULL:
            return CARQUET_OK;

        case CARQUET_LOGICAL_MAP:
        case CARQUET_LOGICAL_LIST:
        case CARQUET_LOGICAL_VARIANT:
            return CARQUET_ERROR_INVALID_ARGUMENT;

        default:
            return CARQUET_ERROR_INVALID_ARGUMENT;
    }
}

carquet_status_t carquet_schema_add_column(
    carquet_schema_t* schema,
    const char* name,
    carquet_physical_type_t physical_type,
    const carquet_logical_type_t* logical_type,
    carquet_field_repetition_t repetition,
    int32_t type_length,
    int32_t parent_index) {

    /* schema and name are nonnull per API contract */

    carquet_status_t status = validate_column_logical_type(
        physical_type, logical_type, type_length);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Validate parent_index: -1 or 0 means root, otherwise must be a valid group */
    if (parent_index == -1) {
        parent_index = 0;
    }
    if (parent_index < 0 || parent_index >= schema->num_elements) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }
    /* Parent must be root (index 0) or a group (no physical type) */
    if (parent_index != 0 && schema->elements[parent_index].has_type) {
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    /* Ensure capacity for new element */
    status = schema_ensure_capacity(schema, schema->num_elements + 1);
    if (status != CARQUET_OK) {
        return status;
    }

    /* Add element to schema */
    int32_t elem_idx = schema->num_elements;
    parquet_schema_element_t* elem = &schema->elements[elem_idx];
    memset(elem, 0, sizeof(*elem));

    elem->name = carquet_arena_strdup(&schema->arena, name);
    elem->has_type = true;
    elem->type = physical_type;
    elem->has_repetition = true;
    elem->repetition_type = repetition;
    elem->type_length = type_length;

    if (logical_type) {
        elem->has_logical_type = true;
        elem->logical_type = *logical_type;
    }

    schema->num_elements++;
    schema->parent_indices[elem_idx] = parent_index;
    schema->elements[parent_index].num_children++;

    /* Compute definition and repetition levels by walking the parent chain */
    int16_t def_level = 0;
    int16_t rep_level = 0;

    if (repetition == CARQUET_REPETITION_OPTIONAL) {
        def_level++;
    } else if (repetition == CARQUET_REPETITION_REPEATED) {
        def_level++;
        rep_level++;
    }

    int32_t ancestor = parent_index;
    while (ancestor > 0) {
        carquet_field_repetition_t ancestor_rep = schema->elements[ancestor].repetition_type;
        if (ancestor_rep == CARQUET_REPETITION_OPTIONAL) {
            def_level++;
        } else if (ancestor_rep == CARQUET_REPETITION_REPEATED) {
            def_level++;
            rep_level++;
        }
        ancestor = schema->parent_indices[ancestor];
    }

    /* Track as leaf */
    schema->leaf_indices[schema->num_leaves] = elem_idx;
    schema->max_def_levels[schema->num_leaves] = def_level;
    schema->max_rep_levels[schema->num_leaves] = rep_level;
    schema->num_leaves++;

    return CARQUET_OK;
}

int32_t carquet_schema_add_group(
    carquet_schema_t* schema,
    const char* name,
    carquet_field_repetition_t repetition,
    int32_t parent_index) {

    /* schema and name are nonnull per API contract */
    if (parent_index == -1) {
        parent_index = 0;
    }
    if (parent_index < 0 || parent_index >= schema->num_elements) {
        return -1;
    }
    /* Parent must be root (index 0) or a group (no physical type) */
    if (parent_index != 0 && schema->elements[parent_index].has_type) {
        return -1;
    }

    /* Ensure capacity for new element */
    if (schema_ensure_capacity(schema, schema->num_elements + 1) != CARQUET_OK) {
        return -1;
    }

    int32_t elem_idx = schema->num_elements;
    parquet_schema_element_t* elem = &schema->elements[elem_idx];
    memset(elem, 0, sizeof(*elem));

    elem->name = carquet_arena_strdup(&schema->arena, name);
    elem->has_type = false;  /* Groups don't have a type */
    elem->has_repetition = true;
    elem->repetition_type = repetition;
    elem->num_children = 0;

    schema->num_elements++;
    schema->parent_indices[elem_idx] = parent_index;
    schema->elements[parent_index].num_children++;

    return elem_idx;
}

/* ============================================================================
 * Schema Queries
 * ============================================================================
 */

int32_t carquet_schema_find_column(
    const carquet_schema_t* schema,
    const char* name) {

    /* schema and name are nonnull per API contract */
    /* Simple linear search */
    for (int32_t i = 0; i < schema->num_leaves; i++) {
        int32_t elem_idx = schema->leaf_indices[i];
        if (schema->elements[elem_idx].name &&
            strcmp(schema->elements[elem_idx].name, name) == 0) {
            return i;
        }
    }

    return -1;
}

int32_t carquet_schema_add_variant(
    carquet_schema_t* schema,
    const char* name,
    carquet_field_repetition_t variant_repetition,
    int32_t parent_index) {

    int32_t outer = carquet_schema_add_group(schema, name, variant_repetition, parent_index);
    if (outer < 0) return -1;

    schema->elements[outer].has_logical_type = true;
    schema->elements[outer].logical_type.id = CARQUET_LOGICAL_VARIANT;
    schema->elements[outer].logical_type.params.variant.specification_version = 1;

    carquet_status_t status = carquet_schema_add_column(
        schema, "metadata", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 0, outer);
    if (status != CARQUET_OK) return -1;

    status = carquet_schema_add_column(
        schema, "value", CARQUET_PHYSICAL_BYTE_ARRAY, NULL,
        CARQUET_REPETITION_REQUIRED, 0, outer);
    if (status != CARQUET_OK) return -1;

    return outer;
}

carquet_status_t carquet_schema_set_field_metadata(
    carquet_schema_t* schema,
    int32_t element_index,
    const char* key,
    const char* value) {

    /* schema and key are nonnull per API contract */
    if (element_index <= 0 || element_index >= schema->num_elements) {
        /* Index 0 is the root group; field metadata attaches to real fields. */
        return CARQUET_ERROR_INVALID_ARGUMENT;
    }

    parquet_schema_element_t* elem = &schema->elements[element_index];

    /* Replace on matching key, else append. Strings live in the schema arena. */
    char* key_copy = carquet_arena_strdup(&schema->arena, key);
    char* val_copy = value ? carquet_arena_strdup(&schema->arena, value) : NULL;
    if (!key_copy || (value && !val_copy)) return CARQUET_ERROR_OUT_OF_MEMORY;

    for (int32_t i = 0; i < elem->num_field_metadata; i++) {
        if (elem->field_metadata[i].key &&
            strcmp(elem->field_metadata[i].key, key) == 0) {
            elem->field_metadata[i].value = val_copy;
            return CARQUET_OK;
        }
    }

    parquet_key_value_t* grown = carquet_arena_alloc(
        &schema->arena,
        (size_t)(elem->num_field_metadata + 1) * sizeof(parquet_key_value_t));
    if (!grown) return CARQUET_ERROR_OUT_OF_MEMORY;
    if (elem->num_field_metadata > 0) {
        memcpy(grown, elem->field_metadata,
               (size_t)elem->num_field_metadata * sizeof(parquet_key_value_t));
    }
    grown[elem->num_field_metadata].key = key_copy;
    grown[elem->num_field_metadata].value = val_copy;
    elem->field_metadata = grown;
    elem->num_field_metadata++;
    return CARQUET_OK;
}

int32_t carquet_schema_num_columns(const carquet_schema_t* schema) {
    /* schema is nonnull per API contract */
    return schema->num_leaves;
}

int32_t carquet_schema_num_elements(const carquet_schema_t* schema) {
    /* schema is nonnull per API contract */
    return schema->num_elements;
}

const carquet_schema_node_t* carquet_schema_get_element(
    const carquet_schema_t* schema,
    int32_t index) {

    /* schema is nonnull per API contract */
    if (index < 0 || index >= schema->num_elements) {
        return NULL;
    }

    /* Return pointer to element (cast as schema_node) */
    return (const carquet_schema_node_t*)&schema->elements[index];
}

/* ============================================================================
 * Schema Node Accessors
 * ============================================================================
 */

const char* carquet_schema_node_name(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return elem->name ? elem->name : "";
}

bool carquet_schema_node_is_leaf(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return elem->has_type;
}

carquet_physical_type_t carquet_schema_node_physical_type(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return elem->type;
}

const carquet_logical_type_t* carquet_schema_node_logical_type(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return elem->has_logical_type ? &elem->logical_type : NULL;
}

carquet_field_repetition_t carquet_schema_node_repetition(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return elem->repetition_type;
}

int16_t carquet_schema_node_max_def_level(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    /* This returns only this node's direct contribution.
     * For accumulated levels, use carquet_schema_max_def_level(). */
    if (elem->repetition_type == CARQUET_REPETITION_OPTIONAL ||
        elem->repetition_type == CARQUET_REPETITION_REPEATED) {
        return 1;
    }
    return 0;
}

int16_t carquet_schema_node_max_rep_level(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return (elem->repetition_type == CARQUET_REPETITION_REPEATED) ? 1 : 0;
}

int32_t carquet_schema_node_type_length(const carquet_schema_node_t* node) {
    /* node is nonnull per API contract */
    const parquet_schema_element_t* elem = (const parquet_schema_element_t*)node;
    return elem->type_length;
}

/* ============================================================================
 * Schema-Level Accessors (accumulated levels for leaf columns)
 * ============================================================================
 */

int16_t carquet_schema_max_def_level(
    const carquet_schema_t* schema,
    int32_t leaf_index) {

    /* schema is nonnull per API contract */
    if (leaf_index < 0 || leaf_index >= schema->num_leaves) {
        return -1;
    }
    return schema->max_def_levels[leaf_index];
}

int16_t carquet_schema_max_rep_level(
    const carquet_schema_t* schema,
    int32_t leaf_index) {

    /* schema is nonnull per API contract */
    if (leaf_index < 0 || leaf_index >= schema->num_leaves) {
        return -1;
    }
    return schema->max_rep_levels[leaf_index];
}

const char* carquet_schema_column_name(
    const carquet_schema_t* schema,
    int32_t leaf_index) {

    /* schema is nonnull per API contract */
    if (leaf_index < 0 || leaf_index >= schema->num_leaves) {
        return NULL;
    }
    int32_t elem_idx = schema->leaf_indices[leaf_index];
    return schema->elements[elem_idx].name;
}

carquet_physical_type_t carquet_schema_column_type(
    const carquet_schema_t* schema,
    int32_t leaf_index) {

    /* schema is nonnull per API contract */
    if (leaf_index < 0 || leaf_index >= schema->num_leaves) {
        return CARQUET_PHYSICAL_BOOLEAN; /* safe default */
    }
    int32_t elem_idx = schema->leaf_indices[leaf_index];
    return schema->elements[elem_idx].type;
}

int32_t carquet_schema_column_path(
    const carquet_schema_t* schema,
    int32_t leaf_index,
    const char** path_out,
    int32_t max_depth) {

    /* schema and path_out are nonnull per API contract */
    if (leaf_index < 0 || leaf_index >= schema->num_leaves || max_depth <= 0) {
        return 0;
    }

    /* Walk from leaf to root, collecting names (excluding root "schema") */
    const char* components[64];
    int32_t depth = 0;

    int32_t elem_idx = schema->leaf_indices[leaf_index];
    while (elem_idx > 0 && depth < 64) {
        components[depth++] = schema->elements[elem_idx].name;
        elem_idx = schema->parent_indices[elem_idx];
    }

    /* Reverse into output (root-first order) */
    int32_t result_len = depth < max_depth ? depth : max_depth;
    for (int32_t i = 0; i < result_len; i++) {
        path_out[i] = components[depth - 1 - i];
    }

    return result_len;
}

/* ============================================================================
 * LIST / MAP Schema Helpers
 * ============================================================================
 */

int32_t carquet_schema_add_list(
    carquet_schema_t* schema,
    const char* name,
    carquet_physical_type_t element_type,
    const carquet_logical_type_t* element_logical_type,
    carquet_field_repetition_t list_repetition,
    int32_t type_length,
    int32_t parent_index) {

    /* Create the outer group with LIST annotation:
     *   <name> (<list_repetition>, LIST) {
     *     list (REPEATED) {
     *       element (<element_type>)
     *     }
     *   }
     */

    /* Outer group: the list container */
    int32_t outer = carquet_schema_add_group(schema, name, list_repetition, parent_index);
    if (outer < 0) return -1;

    /* Set LIST logical type on the outer group */
    schema->elements[outer].has_logical_type = true;
    schema->elements[outer].logical_type.id = CARQUET_LOGICAL_LIST;

    /* Inner repeated group "list" */
    int32_t inner = carquet_schema_add_group(schema, "list", CARQUET_REPETITION_REPEATED, outer);
    if (inner < 0) return -1;

    /* Element leaf column */
    carquet_status_t status = carquet_schema_add_column(
        schema, "element", element_type, element_logical_type,
        CARQUET_REPETITION_OPTIONAL, type_length, inner);
    if (status != CARQUET_OK) return -1;

    return outer;
}

int32_t carquet_schema_add_map(
    carquet_schema_t* schema,
    const char* name,
    carquet_physical_type_t key_type,
    const carquet_logical_type_t* key_logical_type,
    int32_t key_type_length,
    carquet_physical_type_t value_type,
    const carquet_logical_type_t* value_logical_type,
    int32_t value_type_length,
    carquet_field_repetition_t map_repetition,
    int32_t parent_index) {

    /* Create the standard MAP schema:
     *   <name> (<map_repetition>, MAP) {
     *     key_value (REPEATED) {
     *       key (REQUIRED, <key_type>)
     *       value (OPTIONAL, <value_type>)
     *     }
     *   }
     */

    /* Outer group: the map container */
    int32_t outer = carquet_schema_add_group(schema, name, map_repetition, parent_index);
    if (outer < 0) return -1;

    /* Set MAP logical type on the outer group */
    schema->elements[outer].has_logical_type = true;
    schema->elements[outer].logical_type.id = CARQUET_LOGICAL_MAP;

    /* Inner repeated group "key_value" */
    int32_t kv = carquet_schema_add_group(schema, "key_value", CARQUET_REPETITION_REPEATED, outer);
    if (kv < 0) return -1;

    /* Key column (always required) */
    carquet_status_t status = carquet_schema_add_column(
        schema, "key", key_type, key_logical_type,
        CARQUET_REPETITION_REQUIRED, key_type_length, kv);
    if (status != CARQUET_OK) return -1;

    /* Value column (optional) */
    status = carquet_schema_add_column(
        schema, "value", value_type, value_logical_type,
        CARQUET_REPETITION_OPTIONAL, value_type_length, kv);
    if (status != CARQUET_OK) return -1;

    return outer;
}

int32_t carquet_schema_add_list_group(
    carquet_schema_t* schema,
    const char* name,
    carquet_field_repetition_t list_repetition,
    int32_t parent_index) {

    /* Outer LIST-annotated container. */
    int32_t outer = carquet_schema_add_group(schema, name, list_repetition, parent_index);
    if (outer < 0) return -1;
    schema->elements[outer].has_logical_type = true;
    schema->elements[outer].logical_type.id = CARQUET_LOGICAL_LIST;

    /* Inner REPEATED "list" group; caller adds the single element child. */
    return carquet_schema_add_group(schema, "list", CARQUET_REPETITION_REPEATED, outer);
}

int32_t carquet_schema_add_map_group(
    carquet_schema_t* schema,
    const char* name,
    carquet_field_repetition_t map_repetition,
    int32_t parent_index) {

    /* Outer MAP-annotated container. */
    int32_t outer = carquet_schema_add_group(schema, name, map_repetition, parent_index);
    if (outer < 0) return -1;
    schema->elements[outer].has_logical_type = true;
    schema->elements[outer].logical_type.id = CARQUET_LOGICAL_MAP;

    /* Inner REPEATED "key_value" group; caller adds key (required) + value. */
    return carquet_schema_add_group(schema, "key_value", CARQUET_REPETITION_REPEATED, outer);
}

/* ============================================================================
 * Nested Data Helpers
 * ============================================================================
 */

int64_t carquet_count_rows(
    const int16_t* rep_levels,
    int64_t num_values) {

    if (!rep_levels || num_values <= 0) {
        return num_values > 0 ? num_values : 0;
    }

    int64_t rows = 0;
    for (int64_t i = 0; i < num_values; i++) {
        if (rep_levels[i] == 0) rows++;
    }
    return rows;
}

int64_t carquet_list_offsets(
    const int16_t* rep_levels,
    int64_t num_values,
    int16_t list_rep_level,
    int64_t* offsets_out,
    int64_t max_offsets) {

    /* rep_levels and offsets_out are nonnull per API contract */
    if (num_values <= 0 || max_offsets <= 0) {
        return 0;
    }

    /* offsets_out is an Arrow-style offsets array:
     * offsets[i] = start index of list element i
     * offsets[num_lists] = num_values (one past the last)
     * Number of lists = num entries where rep_level < list_rep_level */
    int64_t num_lists = 0;
    for (int64_t i = 0; i < num_values; i++) {
        if (rep_levels[i] < list_rep_level) {
            if (num_lists < max_offsets) {
                offsets_out[num_lists] = i;
            }
            num_lists++;
        }
    }

    /* Write the final offset (one past end) */
    if (num_lists < max_offsets) {
        offsets_out[num_lists] = num_values;
    }

    return num_lists;
}
