# MatBYIB: A MATLAB-Based Code for Bayesian Inference of Inspiral Binary Waveforms with Arbitrary Eccentricity

## Overview
MatBYIB is a MATLAB-based software package designed for Bayesian inference of gravitational wave parameters from inspiraling binary systems with arbitrary eccentricity. This tool is particularly useful for researchers working in gravitational wave astronomy, especially those interested in analyzing extreme mass-ratio inspiral (EMRI) systems using MATLAB. MatBYIB leverages the analytical kludge (AK) waveform model to generate gravitational waveforms and employs the Metropolis-Hastings (M-H) algorithm to perform Bayesian parameter estimation.

## Features
- **Waveform Generation**: Utilizes the AK waveform model to generate gravitational waveforms for EMRI systems with arbitrary eccentricity.
- **Parameter Estimation**: Implements Bayesian inference using the M-H algorithm to estimate parameters such as masses, eccentricity, and spin.
- **Convergence Diagnostics**: Incorporates the Gelman-Rubin convergence criterion to ensure the reliability of the MCMC sampling process.
- **User-Friendly**: Designed with a simple and intuitive MATLAB interface, making it accessible for users familiar with MATLAB.
- **Efficient Sampling**: Employs parallelized MCMC chains to enhance computational efficiency and robustness.
- **Visualization**: Provides tools for visualizing waveforms, posterior distributions, and convergence diagnostics.

## Installation
1. **Clone the Repository**:
   ```bash
   git clone https://github.com/GenliangLi/MatBYIB.git
   ```
2. **Add to MATLAB Path**:
   - Open MATLAB.
   - Navigate to the cloned `MatBYIB` directory.
   - Add the directory to the MATLAB path using the command:
     ```matlab
     addpath('path_to_MatBYIB');
     ```

## Usage

### 1. **Input Parameters**
MatBYIB uses an `input.txt` file to read the parameters of the binary system. The file should contain the following parameters, each on a new line in the format `parameter_name = value`.

#### **Required Input Parameters**
- **`m1`**: Mass of the central black hole (in solar masses, \(M_\odot\)).
  - **Example**: `m1 = 1.0e6` (for \(10^6 M_\odot\))

- **`m2`**: Mass of the orbiting object (in solar masses, \(M_\odot\)).
  - **Example**: `m2 = 10` (for \(10 M_\odot\))

- **`eLSO`**: Eccentricity at the last stable orbit.
  - **Example**: `eLSO = 0.3`

- **`S`**: Spin parameter of the central black hole (\(S/M^2\)).
  - **Example**: `S = 0.01`

- **`z`**: Redshift.
  - **Example**: `z = 0.01`

- **`lambda`**: Angle between the orbital angular momentum and the spin of the central black hole (\(\lambda = \cos^{-1}(\hat{L} \cdot \hat{S})\)).
  - **Example**: `lambda = 60` (in degrees)

- **`Phi_LSO`**: Mean anomaly at the last stable orbit.
  - **Example**: `Phi_LSO = 0` (in degrees)

- **`Gamma_LSO`**: Angle between \(\hat{L} \times \hat{S}\) and the pericenter in the orbital plane.
  - **Example**: `Gamma_LSO = 60` (in degrees)

- **`alpha_LSO`**: Azimuthal direction of \(\hat{L}\) in the orbital plane.
  - **Example**: `alpha_LSO = 60` (in degrees)

- **`thetaS`**: Polar angle of the source direction.
  - **Example**: `thetaS = 60` (in degrees)

- **`phiS`**: Azimuthal angle of the source direction.
  - **Example**: `phiS = 60` (in degrees)

- **`thetaK`**: Polar angle of the central black hole's spin.
  - **Example**: `thetaK = 60` (in degrees)

- **`phiK`**: Azimuthal angle of the central black hole's spin.
  - **Example**: `phiK = 60` (in degrees)

- **`phi0`**: Initial phase of the waveform.
  - **Example**: `phi0 = 0` (in radians)

- **`t_max`**: Total evolution time.
  - **Example**: `t_max = 3.14e6` (in seconds)

- **`detector`**: Detector selection (e.g., LISA, Taiji, Tianqin).
  - **Example**: `detector = LISA`

### 2. **Run the Analysis**
1. **Prepare the Input File**:
   - Save the parameters in the `input.txt` file as shown above.

2. **Run the Main Script**:
   - Open MATLAB and navigate to the MatBYIB directory.
   - Run the `main.m` script:
     ```matlab
     run('main.m');
     ```

3. **View the Results**:
   - The software will generate and display the gravitational waveforms, posterior distributions, and convergence diagnostics.
   - Results will be saved in the specified output directory.

### 3. **Output**
- **Waveform Plots**:
  - Time-domain and frequency-domain gravitational waveforms.
- **Posterior Distributions**:
  - Posterior distributions of the parameters (e.g., masses, eccentricity, spin).
- **Convergence Diagnostics**:
  - Gelman-Rubin convergence factors (\(\hat{R}\)) for each parameter.

### 4. **Further Customization**
- **Adjust MCMC Settings**:
  - Modify the MCMC settings in `MCMC_settings.m` to adjust the number of iterations, initial covariance matrix, and other parameters.
- **Use Different Waveform Models**:
  - The current implementation uses the AK waveform model. For higher accuracy, consider integrating other waveform models (e.g., FastEMRI).

## References
- Genliang Li, Shujie Zhao, Huaike Guo, Zhenheng Li. "MatBYIB: A Matlab-based code for Bayesian inference of inspiral binary waveforms with arbitrary eccentricity." *Journal Not Specified* 2025. [doi: 10.3390/1010000](https://doi.org/10.3390/1010000)


## Contact
For any questions or issues, please contact the authors at:
- Huaike Guo: guohuaike@ucas.ac.cn
- Genliang Li: 2932455541@qq.com

## Documentation
For more detailed information on the software architecture, theory, and numerical examples, please refer to the [MatBYIB GitHub repository](https://github.com/GenliangLi/MatBYIB) and the accompanying [publication](https://doi.org/10.3390/1010000).

