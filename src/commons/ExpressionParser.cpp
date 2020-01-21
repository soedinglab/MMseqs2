#include "ExpressionParser.h"
#include <set>
#include <cstddef>

ExpressionParser::ExpressionParser(const char* expression) : ExpressionParser(expression, {}) {
}

ExpressionParser::ExpressionParser(const char* expression, const std::vector<te_variable>& lookup) {
#define str2(s) #s
#define str(s) str2(s)
#define V(x) { "$" str(x), &variables[((x) - 1)], TE_VARIABLE, NULL},
    vars = {
            V(1) V(2) V(3) V(4) V(5) V(6) V(7) V(8) V(9) V(10) V(11) V(12) V(13) V(14) V(15) V(16)
            V(17) V(18) V(19) V(20) V(21) V(22) V(23) V(24) V(25) V(26) V(27) V(28) V(29) V(30) V(31) V(32)
            V(33) V(34) V(35) V(36) V(37) V(38) V(39) V(40) V(41) V(42) V(43) V(44) V(45) V(46) V(47) V(48)
            V(49) V(50) V(51) V(52) V(53) V(54) V(55) V(56) V(57) V(58) V(59) V(60) V(61) V(62) V(63) V(64)
            V(65) V(66) V(67) V(68) V(69) V(70) V(71) V(72) V(73) V(74) V(75) V(76) V(77) V(78) V(79) V(80)
            V(81) V(82) V(83) V(84) V(85) V(86) V(87) V(88) V(89) V(90) V(91) V(92) V(93) V(94) V(95) V(96)
            V(97) V(98) V(99) V(100) V(101) V(102) V(103) V(104) V(105) V(106) V(107) V(108) V(109) V(110) V(111) V(112)
            V(113) V(114) V(115) V(116) V(117) V(118) V(119) V(120) V(121) V(122) V(123) V(124) V(125) V(126) V(127) V(128)
    };
#undef V
#undef str
#undef str2
    vars.insert(vars.begin(), lookup.begin(), lookup.end());
    expr = te_compile(expression, vars.data(), vars.size(), &err);
}

std::vector<int> ExpressionParser::findBindableIndices() {
    std::vector<int> indices;
    double* base = variables + 0;
    std::vector<const double*> pointers;
    findBound(expr, 0, pointers);
    for (size_t i = 0; i < pointers.size(); ++i) {
        indices.emplace_back(pointers[i] - base);
    }
    std::set<int> s(indices.begin(), indices.end());
    indices.assign(s.begin(), s.end());
    return indices;
}

void ExpressionParser::findBound(const te_expr *n, int depth, std::vector<const double*> &bound) {
#define TYPE_MASK(TYPE) ((TYPE)&0x0000001F)
#define ARITY(TYPE) ( ((TYPE) & (TE_FUNCTION0 | TE_CLOSURE0)) ? ((TYPE) & 0x00000007) : 0 )
    switch(TYPE_MASK(n->type)) {
        case 1:
            break;
        case 0:
            bound.emplace_back(n->bound);
            break;

        default:
            int arity = ARITY(n->type);
            for (int i = 0; i < arity; i++) {
                findBound((const te_expr *) n->parameters[i], depth + 1, bound);
            }
            break;
    }
#undef TYPE_MASK
#undef ARITY
}
