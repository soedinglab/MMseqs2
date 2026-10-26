/**
 * @file test_arrow_nested.c
 * @brief Tests for arbitrary-depth nested Arrow C Data interchange.
 *
 * Exercises the generic Dremel shredder (carquet_writer_write_arrow) and the
 * reassembler (carquet_reader_read_arrow) over struct / list / map composed to
 * several levels, with byte-exact assertions on the reassembled Arrow buffers,
 * nested schema import, error/rejection paths, and page-filter composition with
 * a repeated column.
 */

#include <carquet/carquet.h>
#include "test_helpers.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* ---- minimal Arrow C struct builders (fuzzer-owned, freed by release) ---- */
static void rel_schema(struct ArrowSchema* s) {
    if (!s || !s->release) return;
    free((void*)s->format); free((void*)s->name);
    for (int64_t i = 0; i < s->n_children; i++)
        if (s->children[i]) { if (s->children[i]->release) s->children[i]->release(s->children[i]); free(s->children[i]); }
    free(s->children); s->release = NULL;
}
static void rel_array(struct ArrowArray* a) {
    if (!a || !a->release) return;
    if (a->buffers) { for (int64_t i = 0; i < a->n_buffers; i++) free((void*)a->buffers[i]); free(a->buffers); }
    for (int64_t i = 0; i < a->n_children; i++)
        if (a->children[i]) { if (a->children[i]->release) a->children[i]->release(a->children[i]); free(a->children[i]); }
    free(a->children); a->release = NULL;
}
static char* dups(const char* s){ char* o=malloc(strlen(s)+1); strcpy(o,s); return o; }
static struct ArrowSchema* S(const char* fmt, const char* name, int nullable, int nch){
    struct ArrowSchema* s=calloc(1,sizeof(*s));
    s->format=dups(fmt); s->name=dups(name); s->flags=nullable?ARROW_FLAG_NULLABLE:0;
    s->n_children=nch; s->children=nch?calloc(nch,sizeof(void*)):NULL; s->release=rel_schema; return s;
}
static struct ArrowArray* A(int64_t len,int nbuf,int nch){
    struct ArrowArray* a=calloc(1,sizeof(*a));
    a->length=len; a->null_count=-1; a->n_buffers=nbuf; a->n_children=nch;
    a->buffers=nbuf?calloc(nbuf,sizeof(void*)):NULL; a->children=nch?calloc(nch,sizeof(void*)):NULL;
    a->release=rel_array; return a;
}
static int32_t* i32(const int32_t* x,int64_t n){ int32_t* o=malloc((n?n:1)*4); memcpy(o,x,n*4); return o; }

/* ---- Test: list<list<int32>> round-trip with byte-exact buffer checks ---- */
static int test_list_list(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;

    /* schema: struct { matrix: list<list<int32>> } (nullable lists, req ints) */
    struct ArrowSchema* sc = S("+s","schema",0,1);
    sc->children[0]=S("+l","matrix",1,1);
    sc->children[0]->children[0]=S("+l","element",1,1);
    sc->children[0]->children[0]->children[0]=S("i","element",0,0);
    carquet_schema_t* cs=NULL;
    if (carquet_arrow_import_schema(sc,&cs,&err)!=CARQUET_OK) TEST_FAIL("list_list","import");
    free(sc);  /* import released internals; free the heap node we allocated */

    /* data: [[[1,2],[3]], [], [[4]]] */
    int32_t oo[4]={0,2,2,3}, io[4]={0,2,3,4}, iv[4]={1,2,3,4};
    struct ArrowArray* ii=A(4,2,0); ii->buffers[1]=i32(iv,4);
    struct ArrowArray* il=A(3,2,1); il->buffers[1]=i32(io,4); il->children[0]=ii;
    struct ArrowArray* ol=A(3,2,1); ol->buffers[1]=i32(oo,4); ol->children[0]=il;
    struct ArrowArray* top=A(3,1,1); top->children[0]=ol;

    struct ArrowSchema* sc2=S("+s","schema",0,1);
    sc2->children[0]=S("+l","matrix",1,1);
    sc2->children[0]->children[0]=S("+l","element",1,1);
    sc2->children[0]->children[0]->children[0]=S("i","element",0,0);

    char path[512]; carquet_test_temp_path(path,sizeof(path),"nested_ll");
    carquet_writer_t* w=carquet_writer_create(path,cs,NULL,&err);
    if (!w) TEST_FAIL("list_list","create");
    if (carquet_writer_write_arrow(w,top,sc2,&err)!=CARQUET_OK) TEST_FAIL("list_list", err.message);
    free(top); free(sc2);  /* write_arrow released internals; free heap nodes */
    if (carquet_writer_close(w)!=CARQUET_OK) TEST_FAIL("list_list","close");
    carquet_schema_free(cs);

    /* read back */
    carquet_reader_t* r=carquet_reader_open(path,NULL,&err);
    if (!r) TEST_FAIL("list_list","open");
    struct ArrowSchema as; struct ArrowArray aa;
    if (carquet_reader_read_arrow(r,0,&as,&aa,&err)!=CARQUET_OK) TEST_FAIL("list_list", err.message);

    assert(aa.length==3 && aa.n_children==1);
    struct ArrowArray* outer=aa.children[0];
    assert(strcmp(as.children[0]->format,"+l")==0);
    assert(outer->length==3);
    const int32_t* oofs=(const int32_t*)outer->buffers[1];
    assert(oofs[0]==0 && oofs[1]==2 && oofs[2]==2 && oofs[3]==3);
    struct ArrowArray* inner=outer->children[0];
    assert(inner->length==3);
    const int32_t* iofs=(const int32_t*)inner->buffers[1];
    assert(iofs[0]==0 && iofs[1]==2 && iofs[2]==3 && iofs[3]==4);
    struct ArrowArray* leaf=inner->children[0];
    assert(leaf->length==4);
    const int32_t* vals=(const int32_t*)leaf->buffers[1];
    assert(vals[0]==1 && vals[1]==2 && vals[2]==3 && vals[3]==4);
    assert(leaf->null_count==0);

    aa.release(&aa); as.release(&as);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    TEST_PASS("list_list");
    return 0;
}

/* ---- Test: map<string,int32> with a null value round-trip ---- */
static int test_map(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    struct ArrowSchema* sc=S("+s","schema",0,1);
    sc->children[0]=S("+m","m",1,1);
    sc->children[0]->children[0]=S("+s","entries",0,2);
    sc->children[0]->children[0]->children[0]=S("u","key",0,0);
    sc->children[0]->children[0]->children[1]=S("i","value",1,0);
    carquet_schema_t* cs=NULL;
    if (carquet_arrow_import_schema(sc,&cs,&err)!=CARQUET_OK) TEST_FAIL("map","import");
    free(sc);
    /* leaf count 2 (key, value), key required, value optional */
    assert(carquet_schema_num_columns(cs)==2);
    assert(carquet_schema_max_rep_level(cs,0)==1 && carquet_schema_max_rep_level(cs,1)==1);
    assert(carquet_schema_max_def_level(cs,0)==2);   /* map opt + kv repeated */
    assert(carquet_schema_max_def_level(cs,1)==3);   /* + value opt */

    /* data: {"a":1,"b":NULL}, {}, {"c":3} */
    int32_t oo[4]={0,2,2,3}, ko[4]={0,1,2,3}, vv[3]={1,0,3};
    int vp[3]={1,0,1}; uint8_t* vbits=calloc(1,1); for(int i=0;i<3;i++) if(vp[i]) vbits[0]|=(1u<<i);
    char* kd=malloc(3); memcpy(kd,"abc",3);
    struct ArrowArray* ka=A(3,3,0); ka->buffers[1]=i32(ko,4); ka->buffers[2]=kd;
    struct ArrowArray* va=A(3,2,0); va->buffers[0]=vbits; va->buffers[1]=i32(vv,3); va->null_count=1;
    struct ArrowArray* ea=A(3,1,2); ea->children[0]=ka; ea->children[1]=va;
    struct ArrowArray* ma=A(3,2,1); ma->buffers[1]=i32(oo,4); ma->children[0]=ea;
    struct ArrowArray* top=A(3,1,1); top->children[0]=ma;
    struct ArrowSchema* sc2=S("+s","schema",0,1);
    sc2->children[0]=S("+m","m",1,1);
    sc2->children[0]->children[0]=S("+s","entries",0,2);
    sc2->children[0]->children[0]->children[0]=S("u","key",0,0);
    sc2->children[0]->children[0]->children[1]=S("i","value",1,0);

    char path[512]; carquet_test_temp_path(path,sizeof(path),"nested_map");
    carquet_writer_t* w=carquet_writer_create(path,cs,NULL,&err);
    if (!w) TEST_FAIL("map","create");
    if (carquet_writer_write_arrow(w,top,sc2,&err)!=CARQUET_OK) TEST_FAIL("map", err.message);
    free(top); free(sc2);
    if (carquet_writer_close(w)!=CARQUET_OK) TEST_FAIL("map","close");
    carquet_schema_free(cs);

    carquet_reader_t* r=carquet_reader_open(path,NULL,&err);
    struct ArrowSchema as; struct ArrowArray aa;
    if (carquet_reader_read_arrow(r,0,&as,&aa,&err)!=CARQUET_OK) TEST_FAIL("map", err.message);
    assert(strcmp(as.children[0]->format,"+m")==0);
    struct ArrowArray* m=aa.children[0];
    assert(m->length==3);
    const int32_t* mo=(const int32_t*)m->buffers[1];
    assert(mo[0]==0 && mo[1]==2 && mo[2]==2 && mo[3]==3);
    struct ArrowArray* entries=m->children[0];
    struct ArrowArray* keys=entries->children[0];
    struct ArrowArray* valsA=entries->children[1];
    assert(keys->length==3 && valsA->length==3);
    /* keys "a","b","c" */
    const int32_t* kofs=(const int32_t*)keys->buffers[1];
    const char* kbytes=(const char*)keys->buffers[2];
    assert(kofs[3]==3 && kbytes[0]=='a' && kbytes[1]=='b' && kbytes[2]=='c');
    /* value at index 1 is null */
    assert(valsA->null_count==1);
    const uint8_t* vval=(const uint8_t*)valsA->buffers[0];
    assert((vval[0]&1) && !((vval[0]>>1)&1) && ((vval[0]>>2)&1));

    aa.release(&aa); as.release(&as);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    TEST_PASS("map");
    return 0;
}

/* ---- Test: schema import error/rejection paths ---- */
static void noop_rel(struct ArrowSchema* s){ s->release=NULL; }
static int test_errors(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* out=NULL;

    /* map child that is not a 2-field struct -> INVALID_ARGUMENT */
    struct ArrowSchema kid={0}; kid.format="i"; kid.name="x"; kid.release=noop_rel;
    struct ArrowSchema* mkids[1]={&kid};
    struct ArrowSchema mp={0}; mp.format="+m"; mp.name="m"; mp.n_children=1; mp.children=mkids; mp.release=noop_rel;
    struct ArrowSchema* rk[1]={&mp};
    struct ArrowSchema root={0}; root.format="+s"; root.n_children=1; root.children=rk; root.release=noop_rel;
    assert(carquet_arrow_import_schema(&root,&out,&err)==CARQUET_ERROR_INVALID_ARGUMENT);
    assert(out==NULL);

    /* list with zero children -> INVALID_ARGUMENT */
    struct ArrowSchema lp={0}; lp.format="+l"; lp.name="l"; lp.n_children=0; lp.release=noop_rel;
    struct ArrowSchema* rk2[1]={&lp};
    struct ArrowSchema root2={0}; root2.format="+s"; root2.n_children=1; root2.children=rk2; root2.release=noop_rel;
    assert(carquet_arrow_import_schema(&root2,&out,&err)==CARQUET_ERROR_INVALID_ARGUMENT);

    TEST_PASS("errors");
    return 0;
}

/* ---- Test: page filter composes with a repeated column (RG-level prune) ---- */
static int test_filter_repeated(void) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s=carquet_schema_create(NULL);
    carquet_schema_add_column(s,"id",CARQUET_PHYSICAL_INT64,NULL,CARQUET_REPETITION_REQUIRED,0,0);
    carquet_schema_add_list(s,"vals",CARQUET_PHYSICAL_INT32,NULL,CARQUET_REPETITION_OPTIONAL,0,0);

    carquet_writer_options_t opt; carquet_writer_options_init(&opt);
    opt.write_page_index=true;
    char path[512]; carquet_test_temp_path(path,sizeof(path),"nested_filter");
    carquet_writer_t* w=carquet_writer_create(path,s,&opt,&err);
    if (!w) TEST_FAIL("filter_repeated","create");

    /* RG0: id 1..3 vals [[10],[20],[30]] ; RG1: id 100..102 vals [[40],[50],[60]] */
    for (int g=0; g<2; g++) {
        if (g==1 && carquet_writer_new_row_group(w)!=CARQUET_OK) TEST_FAIL("filter_repeated","new_rg");
        int64_t base = g? 100:1;
        int64_t ids[3]={base,base+1,base+2};
        if (carquet_writer_write_batch(w,0,ids,3,NULL,NULL)!=CARQUET_OK) TEST_FAIL("filter_repeated","id");
        int32_t off[4]={0,1,2,3};
        int32_t vals[3]={(int32_t)(base*10),(int32_t)(base*10+1),(int32_t)(base*10+2)};
        if (carquet_writer_write_list_column(w,1,3,off,NULL,vals,NULL,&err)!=CARQUET_OK) TEST_FAIL("filter_repeated", err.message);
    }
    if (carquet_writer_close(w)!=CARQUET_OK) TEST_FAIL("filter_repeated","close");
    carquet_schema_free(s);

    carquet_reader_t* r=carquet_reader_open(path,NULL,&err);
    if (!r) TEST_FAIL("filter_repeated","open");
    carquet_batch_reader_config_t cfg; carquet_batch_reader_config_init(&cfg);
    carquet_batch_reader_t* br=carquet_batch_reader_create(r,&cfg,&err);
    if (!br) TEST_FAIL("filter_repeated","br");

    /* id >= 100 prunes RG0 entirely. */
    int64_t threshold=100;
    carquet_filter_clause_t clause; memset(&clause,0,sizeof(clause));
    clause.column_index=0; clause.op=CARQUET_FILTER_GE; clause.value=&threshold; clause.value_size=8;
    if (carquet_batch_reader_set_page_filter(br,&clause,1)!=CARQUET_OK) TEST_FAIL("filter_repeated","set_filter");

    carquet_row_batch_t* batch=NULL;
    carquet_status_t st=carquet_batch_reader_next(br,&batch);
    if (st!=CARQUET_OK) TEST_FAIL("filter_repeated","next");
    assert(carquet_row_batch_num_rows(batch)==3);
    const void* iddata=NULL; const uint8_t* nb=NULL; int64_t n=0;
    if (carquet_row_batch_column(batch,0,&iddata,&nb,&n)!=CARQUET_OK) TEST_FAIL("filter_repeated","col");
    const int64_t* ids=(const int64_t*)iddata;
    assert(n==3 && ids[0]==100 && ids[1]==101 && ids[2]==102);  /* RG0 skipped */

    /* No more matching row groups. */
    st=carquet_batch_reader_next(br,&batch);
    assert(st==CARQUET_ERROR_END_OF_DATA);

    carquet_batch_reader_free(br);
    carquet_reader_close(r);
    carquet_test_cleanup(path);
    TEST_PASS("filter_repeated");
    return 0;
}

int main(void) {
    if (test_list_list()) return 1;
    if (test_map()) return 1;
    if (test_errors()) return 1;
    if (test_filter_repeated()) return 1;
    printf("\nAll nested Arrow tests passed.\n");
    return 0;
}
