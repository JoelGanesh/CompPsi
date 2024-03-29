# CompPsi

CompPsi is an application written in C++ that can be used to compute values of $\psi$, where
$$\psi(N) = \sum_{p^k \leq N} \log(p).$$

Three methods are supported:
- Brute-force (i.e., naive);
- Elementary;
- Fourier-theoretic (using _Fast Fourier Transformations_).

The elementary approach is based on a method by [Helfgott & Thompson (2023)](https://link.springer.com/article/10.1007/s40993-022-00408-8). The Fourier-theoretic approach is based on a method by [Hirsch, Kessler & Mendlovic (2023)](https://arxiv.org/abs/2212.09857).
## Dependencies

Before using this application, it is necessary to download the latest version of [Boost](https://www.boost.org/users/download/). Make sure to link the Boost files with the application. It is likely that the steps related to the setup are available at their website.

## Usage

It is possible to set the precision level by updating the definition of ```PRECISION``` in ```Types.h```. Its default is set to ```33```, meaning that computations are done using 33 digits of precision.

After setting up the desired precision level, you can run the application. The intended usage is ```<mode> <initial> <multiplier>```, where ```<mode>``` represents the desired method (```B``` for Brute-force, ```E``` for Elementary or ```F``` for Fourier-theoretic).

## License

This software is licensed under the Boost Software License (Version 1.0):

Permission is hereby granted, free of charge, to any person or organization
obtaining a copy of the software and accompanying documentation covered by
this license (the "Software") to use, reproduce, display, distribute,
execute, and transmit the Software, and to prepare derivative works of the
Software, and to permit third-parties to whom the Software is furnished to
do so, all subject to the following:

The copyright notices in the Software and this entire statement, including
the above license grant, this restriction and the following disclaimer,
must be included in all copies of the Software, in whole or in part, and
all derivative works of the Software, unless such copies or derivative
works are solely in the form of machine-executable object code generated by
a source language processor.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
