corrosion_tests_add_test(parse_target_triple "")
corrosion_tests_add_test(parse_target_triple_should_fail "")

set_tests_properties("parse_target_triple_build" PROPERTIES FAIL_REGULAR_EXPRESSION
        "CMake Warning"
        )

set_tests_properties("parse_target_triple_should_fail_build" PROPERTIES PASS_REGULAR_EXPRESSION
    "CMake Warning"
    )