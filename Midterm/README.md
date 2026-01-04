# Linear ODE General Solution Program

[AI Conversation](https://chatgpt.com/share/695a4d59-b508-800d-aec2-a1610bfdb3f6)

## Introduction
Linear differential equations with constant coefficients appear frequently in mathematics, physics, and engineering, such as in motion models, electrical circuits, and population growth models.
This project presents a Python program that automatically computes the general solution of linear homogeneous ordinary differential equations by transforming the equation into its characteristic polynomial.
The program can handle equations with distinct real roots, repeated real roots, and complex roots.

## Basic Idea
The program solves a differential equation by changing it into a normal polynomial equation.
After converting the equation, the roots of the polynomial are found using the computer.

Different types of roots produce different parts of the general solution.
The program checks the root types and automatically builds the correct solution form.

## How the Program Works
1. The program receives the coefficients of the equation.
2. It forms the characteristic polynomial.
3. It calculates the roots of the polynomial.
4. According to the root types, it constructs the general solution:
  - Real roots → exponential terms
  - Repeated roots → multiply by powers of x
  - Complex roots → sine and cosine terms
5. The program outputs the general solution.

## What Was Learned
Through this project, the relationship between differential equations and polynomial root theory becomes clearer.
It also demonstrates how abstract mathematical formulas can be transformed into practical algorithms that automate symbolic problem solving.

## Conclusion
This program shows that numerical computing tools can be effectively used to solve classical mathematical problems.
By combining root-finding algorithms with theoretical solution rules, the general solution of linear differential equations can be generated quickly and accurately.
