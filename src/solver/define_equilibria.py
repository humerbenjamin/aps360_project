# import standard packages
import numpy as np
import scipy.constants as const
import matplotlib.pyplot as plt

# imports from within this project
from settings import load_machine_geometry, load_solver_parameters
from src.plotting import plot_psi_colormap, plot_psi_colormap_with_contours




def define_multiple_equilibria():

    a=2     # 1.5
    b=0.5   # 0
    c0=1    # 1
    R0=1.5  # 

    # check for non-degeneracy
    if a - c0 <= 0 or c0*R0**2 - b <=0:
        print("degenerate solution")
        return
    else:
        print("non-degenerate solution")

    # calculate the same mesh for all equilibria
    R_1D, Z_1D, R_2D, Z_2D = get_mesh()


    psi_2D = get_psi_on_mesh(R_2D, Z_2D, a, b, c0, R0, show_equilibrium=True)
    psi_lcfs = get_psi_lcfs(psi_2D)
    psi_1D = get_psi_1D(psi_lcfs, 10)
    p_1D = get_p_1D(psi_1D, a)

    plot_psi_colormap_with_contours(R_2D, Z_2D, psi_2D, psi_1D, R0)

    print(psi_1D)

    plt.plot(psi_1D, p_1D)
    plt.savefig("output_plots/test.png")
    plt.close()

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

    return psi_2D


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
def get_psi_1D(psi_lcfs, num_psi):
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
    return np.linspace(0, psi_lcfs, num_psi)


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




####################################################################################################
####################################################################################################
# ASSORTED HELPER FUNCTIONS
####################################################################################################
####################################################################################################

# get the value of psi at the LCFS
####################################################################################################
def get_psi_lcfs(psi_2D, print_lcfs_value=False):
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