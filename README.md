# Calculus C++ Library - c++/discrete math 2025-26

A high-performance calculus library written in C++ with Python bindings. Supports symbolic differentiation and integration for polynomials and trigonometric functions.

## Features

- **Derivatives**: nth-order derivatives for polynomials and trig functions (sin, cos, tan)
- **Integrals**: Indefinite and definite integrals
- **Mixed expressions**: Combine polynomials and trig in one expression
- **Fast evaluation**: C++ backend for performance

## Installation

### Google Colab
```python
!git clone https://github.com/yourusername/your-repo-name.git
%cd your-repo-name
!pip install pybind11
!pip install .
```

### Local Installation
```bash
git clone https://github.com/yourusername/your-repo-name.git
cd your-repo-name
pip install pybind11
pip install .
```

## Usage

See `example.ipynb` for detailed examples.

## API Reference

### Expression Class
- `Expression()` - Create empty expression
- `add_term(coeff, power, trig_type="none")` - Add a term to the expression
  - For polynomials: `add_term(3, 2)` → 3x²
  - For trig: `add_term(2, 0, "sin")` → 2sin(x)
- `evaluate(x)` - Evaluate expression at point x

### Functions
- `derivative(expr, n=1)` - Compute nth derivative (handles all term types)
- `derivative_polynomial(expr, n=1)` - Compute nth derivative of polynomial terms only
- `derivative_trig(expr, n=1)` - Compute nth derivative of trigonometric terms only
- `integral(expr)` - Compute indefinite integral (handles all term types)
- `integral_polynomial(expr)` - Compute indefinite integral of polynomial terms only
- `integral_trig(expr)` - Compute indefinite integral of trigonometric terms only
- `definite_integral(expr, a, b, steps=1000)` - Compute definite integral from a to b using numerical integration

### Supported Functions
- **Polynomials**: x^n (any integer power)
- **Trigonometric**: sin(x), cos(x), tan(x)

## Important Notes

- All trigonometric functions use **radians**, not degrees
- Use `math.radians(degrees)` to convert if needed
- Definite integrals use numerical integration (trapezoidal rule)

## LICENSE

- All of the code is open-source, feel free to use it how you like!
