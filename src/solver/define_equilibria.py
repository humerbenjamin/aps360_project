# import standard packages
    # import packages
import numpy as np
    # import subpackages
import scipy.constants as const
import matplotlib.pyplot as plt
    # import functions
from matplotlib import contour
from skimage import measure


# imports from within this project
from settings import load_machine_geometry, load_solver_parameters
from src.plotting import plot_psi_colormap, plot_psi_colormap_with_contours, plot_P_F_over_psi, plot_n_T_over_psi
from src.load_save import save_single_equilibrium_state


####################################################################################################
####################################################################################################
# MAIN FUNCTION TO DEFINE MULTIPLE EQUILIBRIA
####################################################################################################
####################################################################################################

# main function to define multiple equilibria
####################################################################################################
def define_multiple_equilibria(printouts=False, plot_quantities=False):
    # load in solver parameters
    solver_parameters = load_solver_parameters()
    machine_geometry = load_machine_geometry()

    # initialize linearly spaced arrays
    a_array = np.linspace(solver_parameters['a_min'], solver_parameters['a_max'], solver_parameters['num_a'])
    b_array = np.linspace(solver_parameters['b_min'], solver_parameters['b_max'], solver_parameters['num_b'])
    c0_array = np.linspace(solver_parameters['c0_min'], solver_parameters['c0_max'], solver_parameters['num_c0'])
    R0_array = np.linspace(solver_parameters['R_axis_min'], solver_parameters['R_axis_max'], solver_parameters['num_R_axis'])

    # initialize mesh
    R_1D, Z_1D, R_2D, Z_2D = get_mesh()

    test_params=False
    if test_params:
        a_array=[2]     # 1.5
        b_array=[0.5]   # 0
        c0_array=[1]    # 1
        R0_array=[1.5]  # 

    # calculate states
    statenumber = 0
    degeneracy_counter = 0
    non_degeneracy_counter = 0
    for a in a_array:
        for b in b_array:
            for c0 in c0_array:
                for R0 in R0_array:

                    # check for degeneracy
                    if a - c0 <= 0 or c0*100**2 - b <=0:
                        degeneracy_counter += 1
                        if printouts:
                            print("degenerate solution")
                            print(a, b, c0, R0)

                    # if non degenerate solution 
                    else:
                        if printouts:
                            print("non-degenerate solution")

                        # get 2D psi information
                        psi_2D, psi_axis = get_psi_on_mesh(R_2D, Z_2D, a, b, c0, R0, show_equilibrium=False)

                        if np.min(psi_2D) >= 0:
                            non_degeneracy_counter += 1
                            
                            # additional psi information
                            psi_lcfs, lcfs_contour = find_last_closed_contour(psi_2D)
                            psi_1D = get_psi_1D(psi_lcfs, 0, solver_parameters['num_psi'])
                            psi_contours = get_psi_1D(psi_lcfs, 0, solver_parameters['num_psi_contours'])

                            # other information
                            p_1D = get_p_1D(psi_1D, a)
                            F_1D = get_F_1D(psi_1D, b)
                            n_1D, T_1D = get_n_T_from_p(p_1D, psi_1D)

                            # plot quantities
                            if plot_quantities:
                                plot_psi_colormap_with_contours(R_2D, Z_2D, psi_2D, psi_contours, R0, statenumber)
                                plot_P_F_over_psi(psi_1D, F_1D, p_1D, statenumber)
                                plot_n_T_over_psi(psi_1D, n_1D, T_1D, statenumber)

                            # save data
                            data = {'settings':{'a':a, 'b':b, 'c0':c0, 'R0':R0}, 'geometry':machine_geometry,
                                    'psi_1D':psi_1D, 'F_1D':F_1D, 'p_1D':p_1D, 'n_1D':n_1D, 'T_1D':T_1D
                            }
                            save_single_equilibrium_state(data, statenumber)
                            statenumber += 1
                        
                        else:
                            degeneracy_counter += 1


    print("num_degenerate_solutions:", degeneracy_counter, "// num_non_degenerate_solutions:", non_degeneracy_counter)

    return 


####################################################################################################
####################################################################################################
# DEFINING 2 DIMENSIONAL QUANTITIES
####################################################################################################
####################################################################################################


# get psi on the computed 2-D mesh in R,Z
####################################################################################################
def get_psi_on_mesh(R_2D, Z_2D, a, b, c0, R0, show_equilibrium=False):
    """
    Compute the poloidal flux function (psi) on a 2D R-Z mesh.

    This function evaluates the flux function ψ(R,Z) on a given radial (R) 
    and vertical (Z) meshgrid, based on equilibrium parameters. Optionally, 
    it can also plot the ψ distribution as a colormap for visualization.

    Parameters
    ----------
    R_2D : numpy.ndarray
        2D meshgrid array of radial coordinates (from np.meshgrid).
    Z_2D : numpy.ndarray
        2D meshgrid array of vertical coordinates (from np.meshgrid).
    a : float
        Coefficient controlling the radial dependence of the flux.
    b : float
        Coefficient modifying the vertical dependence of the flux.
    c0 : float
        Coefficient related to shaping of the equilibrium.
    R0 : float
        Reference major radius (center of equilibrium).
    show_equilibrium : bool, optional (default=False)
        If True, produces a colormap plot of ψ(R,Z).

    Returns
    -------
    psi_2D : numpy.ndarray
        2D array of flux values ψ(R,Z), same shape as R_2D and Z_2D.

    (description from ChatGPT)
    """
    psi_2D = np.zeros([len(Z_2D), len(Z_2D[0])])

    for i in range(len(psi_2D)):
        for j in range(len(psi_2D[i])):
            psi_2D[i][j] = psi_of_R_Z(R_2D[i][j], Z_2D[i][j], R0, a, b, c0)

    if show_equilibrium:
        plot_psi_colormap(R_2D, Z_2D, psi_2D)

    psi_axis = psi_of_R_Z(R0, 0, R0, a, b, c0)

    return psi_2D, psi_axis


# compute a 2-D mesh in R,Z 
####################################################################################################
def get_mesh(printout_machine_geometry=False):
    """
    Generate a 2D mesh grid of machine geometry coordinates (R, Z).

    This function loads machine geometry data and constructs a 1D and 2D mesh 
    of the radial (R) and vertical (Z) coordinates. The mesh can be used for 
    plotting or numerical simulations that require spatial grids.

    Parameters
    ----------
    printout_machine_geometry : bool, optional (default=False)
        If True, prints details of the loaded machine geometry for inspection.

    Returns
    -------
    R_1D : numpy.ndarray
        1D array of radial coordinates spanning R_maj_min to R_maj_max.
    Z_1D : numpy.ndarray
        1D array of vertical coordinates spanning Z_min to Z_max.
    R_2D : numpy.ndarray
        2D meshgrid array of radial coordinates (len(Z_1D) x len(R_1D)).
    Z_2D : numpy.ndarray
        2D meshgrid array of vertical coordinates (len(Z_1D) x len(R_1D)).

    (description from ChatGPT)
    """
    machine_geometry_data = load_machine_geometry(printout=printout_machine_geometry)

    R_1D = np.linspace(machine_geometry_data["R_maj_min"],
                       machine_geometry_data["R_maj_max"],
                       machine_geometry_data["num_R_maj"])
    Z_1D = np.linspace(machine_geometry_data["Z_min"],
                       machine_geometry_data["Z_max"],
                       machine_geometry_data["num_Z"])

    R_2D, Z_2D = np.meshgrid(R_1D, Z_1D)
    
    return R_1D, Z_1D, R_2D, Z_2D


####################################################################################################
####################################################################################################
# DEFINING 1 DIMENSIONAL QUANTITIES
####################################################################################################
####################################################################################################


# get 1-D array of psi values from lcfs to the magnetic axis
####################################################################################################
def get_psi_1D(psi_lcfs, psi_axis, num_psi):
    """
    Generate a 1D array of flux values between 0 and the LCFS.

    This function creates a linearly spaced array of ψ values from 
    the magnetic axis (ψ = 0) up to the last closed flux surface (LCFS).

    Parameters
    ----------
    psi_lcfs : float
        Flux value at the last closed flux surface (LCFS).
    num_psi : int
        Number of points in the 1D ψ array.

    Returns
    -------
    psi_1D : numpy.ndarray
        1D array of linearly spaced ψ values from 0 to psi_lcfs.

    (description from ChatGPT)
    """
    return np.linspace(psi_axis, psi_lcfs, num_psi)


# get 1-D array of p values corresponding to psi in a range from the lcfs to the magnetic axis
####################################################################################################
def get_p_1D(psi, a):
    """
    Compute the pressure profile as a function of flux ψ.

    The pressure is defined as a linear function of the distance 
    between each ψ value and the maximum ψ in the array. The 
    expression used is:

        p(ψ) = (a / μ₀) * (ψ_max - ψ)

    where μ₀ is the vacuum permeability.

    Parameters
    ----------
    psi : numpy.ndarray
        1D array of flux values.
    a : float
        Coefficient controlling the pressure scaling.

    Returns
    -------
    p : numpy.ndarray
        1D array of pressure values corresponding to each ψ.

    (description from ChatGPT)
    """
    psimax = np.max(psi)
    p = np.zeros(len(psi))

    for i in range(len(p)):
        p[i] = a / const.mu_0 * (psimax - psi[i])

    return p


# get 1-D array of F values corresponding to psi in a range from the lcfs to the magnetic axis
####################################################################################################
def get_F_1D(psi, b):
    """
    Compute the toroidal field function F(ψ).

    The function is defined as:

        F(ψ) = sqrt(2 * b * ψ)

    where b is a scaling coefficient. This represents the 
    ψ-dependent toroidal field function.

    Parameters
    ----------
    psi : numpy.ndarray
        1D array of flux values.
    b : float
        Coefficient controlling the scaling of F.

    Returns
    -------
    F : numpy.ndarray
        1D array of F(ψ) values corresponding to each ψ.

    (description from ChatGPT)
    """
    psimax = np.max(psi)
    F = np.zeros(len(psi))

    for i in range(len(F)):
        F[i] = np.sqrt(2 * b * psi[i])

    return F


# get 1-D arrays of n and T values corresponding to psi in a range from the lcfs to the magnetic axis
####################################################################################################
def get_n_T_from_p(p_1D, psi_1D, n_plasma=5*10**20):
    n_array = np.full(len(p_1D), n_plasma)
    n_adjusted = np.full(len(p_1D), n_plasma*const.Boltzmann)
    T_array = np.zeros(len(p_1D))
    for i in range(len(T_array)):
        T_array[i] = p_1D[i] / (n_adjusted[i]) * (8.617*10**(-5))

    T_array[-1] = T_array[-2]
    n_array[-1] = 0

    return n_array, T_array


####################################################################################################
####################################################################################################
# ASSORTED HELPER FUNCTIONS
####################################################################################################
####################################################################################################

# get the value of psi at the LCFS
####################################################################################################
def get_psi_lcfs(psi_2D, err=1e-5, num_iter=50, print_lcfs_value=False):
    """
    Determine the flux value at the last closed flux surface (LCFS).

    The LCFS is approximated as the minimum flux value along the 
    outer boundary of the poloidal flux grid ψ(R,Z). This function 
    checks the minimum along the top, bottom, left, and right edges 
    of the 2D array and returns the overall minimum.

    Parameters
    ----------
    psi_2D : numpy.ndarray
        2D array of flux values ψ(R,Z).
    print_lcfs_value : bool, optional (default=False)
        If True, prints the LCFS flux value in LaTeX-style formatting.

    Returns
    -------
    lcfs_psi : float
        The minimum boundary flux value, representing ψ at the LCFS.

    (description from ChatGPT)
    """
    # get minimum on each edge
    top_min = np.min(psi_2D[0, :])
    bot_min = np.min(psi_2D[-1, :])
    left_min = np.min(psi_2D[:, 0])
    right_min = np.min(psi_2D[:, -1])

    # get overall minimum (lcfs)
    lcfs_psi = min(top_min, bot_min, left_min, right_min)
    if print_lcfs_value:
        print(r"$\psi_\mathrm{lcfs}$" + f"={lcfs_psi:.4f}")


    return lcfs_psi


def is_closed(contour, tol=1e-6):
    """Check if a contour is closed (first and last points coincide)."""
    return np.linalg.norm(contour[0] - contour[-1]) < tol


def find_last_closed_contour(arr, start_level=None, end_level=None, step=0.01, verbose=False):
    """
    Iteratively adjusts contour level to find the last closed contour.

    Parameters
    ----------
    arr : 2D numpy array
        Scalar field (e.g. grayscale or elevation map)
    start_level : float, optional
        Starting contour level (defaults to min(arr))
    end_level : float, optional
        Ending contour level (defaults to max(arr))
    step : float
        Increment for level scanning
    verbose : bool
        Print debug info

    Returns
    -------
    level : float
        Level corresponding to the last closed contour
    contour : np.ndarray or None
        Coordinates of the last closed contour
    """
    start_level = start_level if start_level is not None else np.min(arr)
    end_level = end_level if end_level is not None else np.max(arr)

    last_closed = None
    last_level = None

    # Sweep from low → high level
    for level in np.arange(start_level, end_level, step):
        contours = measure.find_contours(arr, level=level)
        closed = [c for c in contours if is_closed(c)]

        if closed:
            last_closed = closed[-1]   # last closed contour at this level
            last_level = level

            if verbose:
                print(f"Level {level:.3f}: {len(closed)} closed contour(s) found")
        else:
            # Once we lose all closed contours, break
            if last_closed is not None:
                break

    return last_level, last_closed



# get the value of psi at a given R, Z pair for constants R0, a, b, c0
####################################################################################################
def psi_of_R_Z(R, Z, R0, a, b, c0):
    """
    Evaluate the poloidal flux function ψ(R, Z).

    The flux ψ is defined as:

        ψ(R, Z) = 0.5 * (c0 * R² - b) * Z²
                  + (a - c0) / 8 * (R² - R0²)²

    Parameters
    ----------
    R : float or numpy.ndarray
        Radial coordinate(s).
    Z : float or numpy.ndarray
        Vertical coordinate(s).
    R0 : float
        Reference major radius.
    a : float
        Coefficient controlling radial dependence.
    b : float
        Coefficient modifying vertical dependence.
    c0 : float
        Coefficient related to shaping of the equilibrium.

    Returns
    -------
    psi : float or numpy.ndarray
        Value(s) of the poloidal flux ψ(R, Z).

    (description from ChatGPT)
    """
    # return c0 * (0.5*(R**2 - R0**2)**2 + 0.5*a**2*Z**2 - 0.25*b*Z**2)
    return 0.5 * (c0 * R**2 - b) * Z**2 + (a - c0) / 8 * (R**2 - R0**2)**2