/* Currently unused subclass of FFTLibrary using the FFTW library.
		class FFTW : public FFTLibrary<double, fftw_complex, fftw_complex*>
		{
			public:
			// Constructor
			FFTW(int size, int n)
			{
				// To compute the non-modular convolution on n vectors,
				// the vectors need to be padded with zeros until
				// each vector has size equal to n times the original size.
				this->size = size * n;

				// Initialize arrays for input and output by allocating memory.
				in = (double*)fftw_malloc(sizeof(double) * this->size);
				out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);

				inv_in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);
				inv_out = (double*)fftw_malloc(sizeof(double) * this->size);

				// Set up plans for both the forward and backward Fourier transform
				// to be executed multiple times.
				forward_plan = fftw_plan_dft_r2c_1d(this->size, in, out, FFTW_ESTIMATE);
				backward_plan = fftw_plan_dft_c2r_1d(this->size, inv_in, inv_out, FFTW_ESTIMATE);
			}

			// Destructor
			~FFTW()
			{
				// Destroy the plans, free the allocated space for in and out, and clean up.
				fftw_destroy_plan(forward_plan);
				fftw_destroy_plan(backward_plan);
				fftw_free(in); fftw_free(out);
				fftw_free(inv_in); fftw_free(inv_out);
				fftw_cleanup();
			}

			private:
			// Serves for real input to compute forward/backward DFTs.
			double* in;
			fftw_complex* inv_in;

			// Serves for complex output to store forward/backward DFTs.
			fftw_complex* out;
			double* inv_out;

			// A plan created to compute FFTs with data from 'in' (real) to 'out' (complex).
			fftw_plan forward_plan;
			fftw_plan backward_plan;

			void InitializeContainer(fftw_complex*& Fv) override
			{
				Fv = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * this->size);
				for (int i = 0; i < size; i++)
				{
					Fv[i][RE] = 1.0;
					Fv[i][IM] = 0.0;
				}
			}

			void DestructContainer(fftw_complex*& Fv) override
			{
				fftw_free(Fv);
			}

			void Multiply(fftw_complex*& Fv, fftw_complex* Fv_i) override
			{
				for (int j = 0; j < size; j++)
				{
					double temp = Fv[j][RE] * Fv_i[j][RE] - Fv[j][IM] * Fv_i[j][IM];

					Fv[j][IM] = Fv[j][RE] * Fv_i[j][IM] + Fv[j][IM] * Fv_i[j][RE];
					Fv[j][RE] = temp;
				}
			}

			// Implementation of FFT using the FFTW library.
			fftw_complex* FFT(std::vector<double> v) override
			{
				for (int j = 0; j < size; j++)
				{
					in[j] = v[j];
				}
				fftw_execute(forward_plan);
				return out;
			}

			// Implementation of IFFT using the FFTW library.
			std::vector<double> IFFT(fftw_complex* Fv) override
			{
				for (int j = 0; j < size; j++)
				{
					inv_in[j][RE] = Fv[j][RE];
					inv_in[j][IM] = Fv[j][IM];
				}
				fftw_execute(backward_plan);

				std::vector<double> result(inv_out, inv_out + size);
				return result;
			}

		};*/