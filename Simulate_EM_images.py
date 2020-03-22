# %%
import mrcfile
from matplotlib import pyplot as plt
import numpy as np
from scipy.constants import physical_constants
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import zoom

# %%
def wavelength(V):
    e = physical_constants['elementary charge'][0] # Coulombs
    m = physical_constants['electron mass'][0] # kg
    h = physical_constants['Planck constant'][0] # J s
    c = physical_constants['speed of light in vacuum'][0] # m/s

    V = V * 1e3 # kV

    # Electron wavelength in Angstroms for 200 keV acceleration potential
    lam = h*c/np.sqrt((e*V)**2 + 2*e*V*m*c**2) * 1e10

    print(f"Accelerating voltage {V} volts, electron wavelength {lam:.2f} Å")
    return lam

# %%
def ctf_2d(box,angpix, defocus_U, defocus_V, astigmatism, amplitude_contrast, lam, spherical_aberration):
    """Compute the CTF.

    Evaluates the CTF at `omega`, returning an array of shape `omega.shape[0:-1]`

    Parameters
    ----------
    omega : `numpy.ndarray` of shape (..., 2) 
        Spatial frequencies at which to evaluate the CTF

        In units of 1/A
    defocus_U, defocus_V : float
        Maximal and minimal defocus
      
        In units of Angstrom
    amplitude_contrast : float
        Amplitude contrast

        Unitless, in [0, 1]
    astigmatism : float
        Astigmatistm angle

        In units of radians
    lam : float
        Electron wavelength

        In units of Angstroms
    spherical_aberration : float
        Spherical aberration

        In units of Angstroms

    Returns
    -------
    `numpy.ndarray` of shape `omega.shape[0:-1]`
        CTF evaluated at requested points

    """
    print('unit check '+str(defocus_U) + ' ' + str(defocus_V)+ ' ' + str(lam) + ' ' + str(spherical_aberration) + '\n')
    
    freqs = np.fft.fftshift(np.fft.fftfreq(int(box), angpix))

    xi_y, xi_x = np.meshgrid(freqs, freqs, indexing='ij')

    omega = np.stack((xi_x, xi_y)).T
    
    omega_abs = np.sqrt(np.sum(omega*omega, axis=-1))

    # Compute defocus term with astigmatism
    thetas = np.arctan(omega[..., 1], omega[..., 0])

    defocus = defocus_U*(np.cos(thetas - astigmatism)**2) + defocus_V*(np.sin(thetas - astigmatism)**2)
        
    gamma =   np.pi*(defocus*lam*(omega_abs**2) - 0.5*(spherical_aberration)*(lam**3)*(omega_abs**4))

    return -np.sqrt(1 - amplitude_contrast**2)*np.sin(gamma) - amplitude_contrast*np.cos(gamma)

# %%
def simulate_image_direct(mrcf, v, V, defocus, C_s):
    """Simulate a pair of images, with and without Ewald Sphere curvature
    effects
    
    Parameters
    ----------
    vol_hat : `scipy.interpolate.RegularGridInterpolator`
        molecule's 3D Fourier Transform
    v : array-like or `numpy.ndarray` of size (3,3)
        image projection direction, in which case a random in-plane
        rotation is assigned, or three basis vectors for the image
        plane
    lam : float
        Electron wavelength in Angstroms
    defocus : float
        Image defocus in Angstroms
    C_s : float
        Spherical aberration in Angstroms
        
    Returns
    -------
    tuple of `numpy.ndarray`
        pair of images, with and without Ewald Sphere curvature effects

    """
    # Read the molecule
    with mrcfile.open(mrcf) as rho_mrc:
        rho_samples = rho_mrc.data

        L = rho_mrc.header.cella.x
        N = rho_mrc.data.shape[0]

    angpix = L/N

    assert N == rho_samples.shape[1] and N == rho_samples.shape[2] 

    print(f"Loaded {mrcf}: {N}^3 voxels, physical dimensions {L:.2f}^3 "
      f"Å^3, {angpix:.2f} Å/voxel")
    if not isinstance(v, np.ndarray) or v.shape == (3,):
        # Generate a rotation matrix from the projection direction
        # v. We'll take a random in-plane rotation.
        a = np.cross(v, 2*np.random.random(3) - 1)
        b = np.cross(v, a)
    
        assert not np.isclose(np.linalg.norm(a), 0), ("We randomly generated "
                                                      "a singular vector "
                                                      "orthogonal to v, "
                                                      "try again")
        a = a / np.linalg.norm(a)
        b = b / np.linalg.norm(b)
        c = v / np.linalg.norm(v)
    
        R = np.array([a,b,c]).T
    else:
        assert v.shape == (3,3), ("Must provide a rotation matrix "
                                  "representing this projection")
        R = v
  
    assert np.isclose(np.linalg.det(R), 1), ("Unexpected: rotation matrix "
                                             "should have determinant 1")
    assert np.isclose(np.dot(R[:,0], R[:,1]), 0), ("Unexpected: a and b "
                                                   "should be perpendicular")
    assert np.isclose(np.dot(R[:,1], R[:,2]), 0), ("Unexpected: b and c "
                                                   "should be perpendicular")
    assert np.isclose(np.dot(R[:,0], R[:,2]), 0), ("Unexpected: a and c "
                                                   "should be perpendicular")
       
        
    #N is box size
    freqs = np.fft.fftshift(np.fft.fftfreq(N, angpix))

    xi_y, xi_x = np.meshgrid(freqs, freqs, indexing='ij')
    
    s2 = xi_y**2 + xi_x**2
    s = np.sqrt(s2)
    
    if N % 2 == 0:
        friedel_1d = np.roll(np.arange(N)[::-1], 1)
    else:
        friedel_1d = np.arange(N)[::-1]
    
    friedel_mate_indices = tuple(np.meshgrid(friedel_1d, friedel_1d, indexing='ij'))

    
    # Calculate defocus term
    lam = wavelength(200) #angstrom
    chi = np.pi*(0.5*C_s*lam**3*s**4 - defocus*lam*s**2)

    freqs_integer = np.arange(-N//2, N//2, dtype=np.int16)
    l, k, h = np.meshgrid(freqs_integer, freqs_integer, freqs_integer, indexing='ij')

    rho_dft = np.fft.fftshift(np.fft.fftn(rho_samples))
    

    phase_flip = (h + k + l) % 2
    phase_flip *= -2
    phase_flip += 1
    rho_hat_samples = angpix**3*phase_flip*rho_dft
    
    vol_hat = RegularGridInterpolator((freqs, freqs, freqs), rho_hat_samples,bounds_error=False, fill_value = 0.0)

    # Also generate 2D integer frequency grids for applying the opposite
    # phaseflip when generating images
    k_2d, h_2d = np.meshgrid(freqs_integer, freqs_integer, indexing='ij')
    
    # Compute image with no Ewald curvature
    xi_rot = np.einsum('ij,jkl->ikl', R,
                       np.stack((xi_x, xi_y, np.zeros((N,N)))))
    
    phase_shift = (h_2d + k_2d) % 2
    phase_shift *= -2
    phase_shift += 1

    I_hat_no_curvature = vol_hat((xi_rot[2,...], xi_rot[1,...], xi_rot[0,...]))
    
    xi_ctf = np.stack((xi_x, xi_y)).T
    #Defocus given as angstrom, but in ctf_2d takes as nm
    ctf_samples = ctf_2d(N,angpix, defocus, defocus, 0, 0.1, lam, C_s) 
    I_hat_no_curvature_ctf = I_hat_no_curvature*ctf_samples

    I_hat_no_curvature = (1/angpix)**2*phase_shift*I_hat_no_curvature
    I_hat_no_curvature_ctf = (1/angpix)**2*phase_shift*I_hat_no_curvature_ctf
    
    # Since we sample the no-curvature image off-axis, it will be
    # interpolated and not perfectly satisfy Freidel symmetry. We
    # average to enforce Freidel symmetry
    I_hat_no_curvature = 0.5*(
        I_hat_no_curvature +
        np.conj(I_hat_no_curvature[friedel_mate_indices]))

    I_hat_no_curvature_ctf = 0.5*(
        I_hat_no_curvature_ctf +
        np.conj(I_hat_no_curvature_ctf[friedel_mate_indices]))
    
    # Correct Nyquist components, see comment above
    #if N % 2 == 0:
    #    I_hat_no_curvature[0,:] = 0 
    #    I_hat_no_curvature[:,0] = 0 
    
    I_no_curvature = np.fft.ifftn(np.fft.ifftshift(I_hat_no_curvature))
    I_no_curvature_ctf = np.fft.ifftn(np.fft.ifftshift(I_hat_no_curvature_ctf))

    imag_magnitude = np.linalg.norm(np.imag(I_no_curvature))
    assert imag_magnitude < 1e-10, ("Image has large "
                                    f"imaginary magnitude {imag_magnitude}")

    imag_magnitude = np.linalg.norm(np.imag(I_no_curvature_ctf))
    assert imag_magnitude < 1e-10, ("Image+CTF has large "
                                    f"imaginary magnitude {imag_magnitude}")

    I_no_curvature = np.real(I_no_curvature)
    I_no_curvature_ctf = np.real(I_no_curvature_ctf)    
    return I_no_curvature, I_no_curvature_ctf, R

# %%
I,I_ctf, R = simulate_image_direct('/Users/ekellogg/projects/datasets/STC/cryosparc_exp000067_006.mrc', [0, 1, 1], 300, 30000, 2.7)

# %%
from matplotlib import pyplot as plt

# %%
#replace with the image you want to test with
#mrcf = mrcfile.open('/Users/ekellogg/projects/datasets/STC/cryosparc_exp000067_006.mrc')
#defocus = 30000
#tt=ctf_2d(mrcf.data.shape[0],1.16,defocus,defocus,0,0.1,wavelength(200),2.7)
#tt.shape

# %%
#plt.imshow(I, cmap='gray')

# %%
#plt.imshow(I_ctf, cmap='gray')

# %%
