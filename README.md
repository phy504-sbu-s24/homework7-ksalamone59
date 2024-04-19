# Homework #7

### due: Fri April 26, 2024

> Each student is responsible for doing their own work.  You are welcome to
> discuss this on slack or in class, but when it comes to writing the code,
> you are expected to write it yourself (and not just copy from someone).

> Note: using ChatGPT or similar AI tools to write your code is not allowed.

1. In our [debugging
   experiment](https://zingale.github.io/phy504/debugging.html), we
   looked at the trapezoid rule for integration.

   Given a function, $f(x)$, defined on $[a, b]$, we divide the domain
   into $N$ slabs, and the integral is approximated as:

   $$I = \frac{1}{2} \Delta x \sum_{i=1}^N \left ( f(x_i) + f(x_{i+1}) \right )$$

   where $\Delta x$ is the width of a slab and $x_i$ and $x_{i+1}$ are the
   coordinates of the left and right edge of slab $i$.

   a. Write a C++ function to integrate a function $f(x)$ using the
      trapezoid rule.  Your C++ function should take as input the number
      of slabs, $N$, and the integrand, $f(x)$, using `std::function`, as
      we saw in class.

   b. Write a second C++ function to do Monte Carlo integration.  The interface
      will be similar to your trapezoid function, but now $N$ will be the number
      of random points.

   c. Write a program that integrates

      $$I = \int_{-5}^5 e^{-x^2} dx$$

      using both the trapezoid rule and Monte Carlo.  Vary N
      (say from 8 to 1024, by doubling), and have your program output
      a table with: $N$, trapezoid answer, MC answer.

   d. Now do the same, but for the integral:

      $$I = \int_0^2 \sin^2 \left ( \frac{1}{x (2-x) + \epsilon} \right ) dx$$

      where $\epsilon \approx 10^{-12}$ is there to prevent dividing by zero.

      You might want to plot this integrand to see the structure.
