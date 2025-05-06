# electronic_component_cooling

# Newton’s Law of Cooling: Analytical and Numerical Approaches

The cooling of electronic components is essential for their proper functioning and to extend their lifespan. This process, in a controlled environment, can be modeled by **Newton’s Law of Cooling**, described by the following differential equation:

$$
\frac{dT}{dt} = -k(T - T_a)
$$

Where:
- \( T(t) \): Temperature of the component at time \( t \) (in °C)
- \( T_a \): Ambient temperature (assumed constant)
- \( k \): Positive constant related to the cooling rate (in \( \text{s}^{-1} \))

## Objective

This project aims to explore the provided differential equation both **analytically** and **numerically**, in order to understand its behavior and the effectiveness of each method.

## Analytical Solution

The analytical solution gives an exact expression for the temperature over time. However, for more complex systems, this form may not be available or practical to obtain. In such cases, **numerical methods** become essential.

## Numerical Methods

Two widely-used numerical techniques will be implemented:

- **Euler’s Method**
- **Fourth-Order Runge-Kutta Method (RK4)**

These methods will be applied using different step sizes to evaluate:
- **Accuracy** in approximating the analytical solution
- **Efficiency** based on computational cost and precision

## Analysis and Comparison

- The results from both numerical methods will be **graphically compared** to the analytical solution.
- **Percentage errors** will be calculated to evaluate each method's precision.
- The influence of **step size** on accuracy will be discussed.

## Conclusion

This study illustrates the application of differential equations to real-world problems, highlighting the importance of numerical methods in **engineering** and **science** when analytical solutions are difficult or impossible to obtain.
