#ifndef EXPRESSION_PARSER_H
#define EXPRESSION_PARSER_H

#include <vector>
#include <tinyexpr.h>

class ExpressionParser {
public:
    ExpressionParser(const char* expression);
    ExpressionParser(const char* expression, const std::vector<te_variable>& lookup);

    ~ExpressionParser() {
        if (expr) {
            te_free(expr);
        }
    }

    bool isOk() {
        return err == 0;
    }

    std::vector<int> findBindableIndices();

    void bind(unsigned int index, double value) {
        if (index > 127) {
            return;
        }
        variables[index] = value;
    }

    double evaluate() {
        return te_eval(expr);
    }

private:
    void findBound(const te_expr *n, int depth, std::vector<const double*> &bound);

    te_expr *expr;
    std::vector<te_variable> vars;
    double variables[128];
    int err;
};

#endif
