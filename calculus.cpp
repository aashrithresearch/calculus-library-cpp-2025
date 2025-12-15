#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>
#include <sstream>

namespace py = pybind11;

// Represents a term in a polynomial or trig expression
struct Term {
    double coeff;
    int power;
    std::string trig_type; // "none", "sin", "cos", "tan"
    
    Term(double c = 0, int p = 0, std::string t = "none") 
        : coeff(c), power(p), trig_type(t) {}
};

class Expression {
private:
    std::vector<Term> terms;
    
public:
    Expression() {}
    Expression(const std::vector<Term>& t) : terms(t) {}
    
    void add_term(double coeff, int power, std::string trig_type = "none") {
        terms.push_back(Term(coeff, power, trig_type));
    }
    
    const std::vector<Term>& get_terms() const { return terms; }
    
    // Evaluate expression at x
    double evaluate(double x) const {
        double result = 0.0;
        for (const auto& term : terms) {
            double val = term.coeff;
            if (term.trig_type == "sin") {
                val *= std::sin(x);
            } else if (term.trig_type == "cos") {
                val *= std::cos(x);
            } else if (term.trig_type == "tan") {
                val *= std::tan(x);
            }
            if (term.power != 0) {
                val *= std::pow(x, term.power);
            }
            result += val;
        }
        return result;
    }
    
    std::string to_string() const {
        if (terms.empty()) return "0";
        std::ostringstream oss;
        bool first = true;
        for (const auto& term : terms) {
            if (!first && term.coeff >= 0) oss << " + ";
            else if (term.coeff < 0) oss << " - ";
            else if (first) oss << "";
            
            double c = std::abs(term.coeff);
            if (c != 1.0 || (term.power == 0 && term.trig_type == "none")) {
                oss << c;
            }
            
            if (term.trig_type != "none") {
                oss << term.trig_type << "(x)";
            }
            
            if (term.power > 0) {
                if (term.trig_type != "none") oss << "*";
                oss << "x";
                if (term.power > 1) oss << "^" << term.power;
            }
            first = false;
        }
        return oss.str();
    }
};

// Polynomial derivative (nth order)
Expression derivative_polynomial(const Expression& expr, int n = 1) {
    if (n < 0) throw std::invalid_argument("Order must be non-negative");
    if (n == 0) return expr;
    
    Expression result;
    for (const auto& term : expr.get_terms()) {
        if (term.trig_type != "none") continue; // Skip trig terms
        
        double coeff = term.coeff;
        int power = term.power;
        
        // Apply power rule n times
        for (int i = 0; i < n; ++i) {
            if (power == 0) {
                coeff = 0;
                break;
            }
            coeff *= power;
            power--;
        }
        
        if (coeff != 0) {
            result.add_term(coeff, power);
        }
    }
    return result;
}

// Trigonometric derivative
Expression derivative_trig(const Expression& expr, int n = 1) {
    if (n < 0) throw std::invalid_argument("Order must be non-negative");
    if (n == 0) return expr;
    
    Expression result;
    for (const auto& term : expr.get_terms()) {
        if (term.trig_type == "none") continue; // Skip polynomial terms
        
        double coeff = term.coeff;
        std::string trig = term.trig_type;
        
        // Apply chain rule n times
        for (int i = 0; i < n; ++i) {
            if (trig == "sin") {
                trig = "cos";
            } else if (trig == "cos") {
                trig = "sin";
                coeff = -coeff;
            } else if (trig == "tan") {
                // d/dx(tan(x)) = sec^2(x) - simplified as 1 + tan^2(x)
                throw std::runtime_error("Higher order tan derivatives not implemented");
            }
        }
        
        result.add_term(coeff, term.power, trig);
    }
    return result;
}

// General derivative (handles both polynomial and trig)
Expression derivative(const Expression& expr, int n = 1) {
    Expression poly_result = derivative_polynomial(expr, n);
    Expression trig_result = derivative_trig(expr, n);
    
    Expression combined;
    for (const auto& term : poly_result.get_terms()) {
        combined.add_term(term.coeff, term.power, term.trig_type);
    }
    for (const auto& term : trig_result.get_terms()) {
        combined.add_term(term.coeff, term.power, term.trig_type);
    }
    return combined;
}

// Polynomial integral
Expression integral_polynomial(const Expression& expr) {
    Expression result;
    for (const auto& term : expr.get_terms()) {
        if (term.trig_type != "none") continue; // Skip trig terms
        
        int new_power = term.power + 1;
        double new_coeff = term.coeff / new_power;
        result.add_term(new_coeff, new_power);
    }
    return result;
}

// Trigonometric integral
Expression integral_trig(const Expression& expr) {
    Expression result;
    for (const auto& term : expr.get_terms()) {
        if (term.trig_type == "none") continue; // Skip polynomial terms
        
        double coeff = term.coeff;
        std::string new_trig;
        
        if (term.trig_type == "sin") {
            new_trig = "cos";
            coeff = -coeff;
        } else if (term.trig_type == "cos") {
            new_trig = "sin";
        } else if (term.trig_type == "tan") {
            throw std::runtime_error("tan integral not implemented in simple form");
        }
        
        result.add_term(coeff, term.power, new_trig);
    }
    return result;
}

// General integral (handles both polynomial and trig)
Expression integral(const Expression& expr) {
    Expression poly_result = integral_polynomial(expr);
    Expression trig_result = integral_trig(expr);
    
    Expression combined;
    for (const auto& term : poly_result.get_terms()) {
        combined.add_term(term.coeff, term.power, term.trig_type);
    }
    for (const auto& term : trig_result.get_terms()) {
        combined.add_term(term.coeff, term.power, term.trig_type);
    }
    return combined;
}

// Definite integral using numerical integration (trapezoidal rule)
double definite_integral(const Expression& expr, double a, double b, int steps = 1000) {
    double h = (b - a) / steps;
    double sum = 0.5 * (expr.evaluate(a) + expr.evaluate(b));
    
    for (int i = 1; i < steps; ++i) {
        sum += expr.evaluate(a + i * h);
    }
    
    return sum * h;
}

PYBIND11_MODULE(calculus, m) {
    m.doc() = "Calculus library with derivatives and integrals";
    
    py::class_<Term>(m, "Term")
        .def(py::init<double, int, std::string>(), 
             py::arg("coeff") = 0, py::arg("power") = 0, py::arg("trig_type") = "none")
        .def_readwrite("coeff", &Term::coeff)
        .def_readwrite("power", &Term::power)
        .def_readwrite("trig_type", &Term::trig_type);
    
    py::class_<Expression>(m, "Expression")
        .def(py::init<>())
        .def(py::init<const std::vector<Term>&>())
        .def("add_term", &Expression::add_term, 
             py::arg("coeff"), py::arg("power"), py::arg("trig_type") = "none")
        .def("evaluate", &Expression::evaluate)
        .def("__str__", &Expression::to_string)
        .def("__repr__", &Expression::to_string);
    
    m.def("derivative", &derivative, 
          py::arg("expr"), py::arg("n") = 1,
          "Compute nth derivative of expression");
    
    m.def("derivative_polynomial", &derivative_polynomial,
          py::arg("expr"), py::arg("n") = 1,
          "Compute nth derivative of polynomial terms");
    
    m.def("derivative_trig", &derivative_trig,
          py::arg("expr"), py::arg("n") = 1,
          "Compute nth derivative of trigonometric terms");
    
    m.def("integral", &integral,
          py::arg("expr"),
          "Compute indefinite integral of expression");
    
    m.def("integral_polynomial", &integral_polynomial,
          py::arg("expr"),
          "Compute indefinite integral of polynomial terms");
    
    m.def("integral_trig", &integral_trig,
          py::arg("expr"),
          "Compute indefinite integral of trigonometric terms");
    
    m.def("definite_integral", &definite_integral,
          py::arg("expr"), py::arg("a"), py::arg("b"), py::arg("steps") = 1000,
          "Compute definite integral using numerical integration");
}
