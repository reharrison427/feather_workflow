import glob
import numpy as np
import os

from astropy.io import fits
import matplotlib.pyplot as plt

'''
Feathers TP and interferometer data and makes diagnostic plots 
showing proportion of TP flux recovered in feathered image.
'''


def casa_to_fits(casa_image, fits_file=None, dropdeg=False):
    """
    Export a CASA image/PB to FITS.

    If fits_file is not supplied:
      something.image -> something.fits
      something.pb    -> something.pb.fits
    """

    if fits_file is None:
        if casa_image.endswith('.image'):
            fits_file = casa_image[:-6] + '.fits'
        else:
            fits_file = casa_image + '.fits'

    exportfits(
        imagename=casa_image,
        fitsimage=fits_file,
        overwrite=True,
        dropdeg=dropdeg
    )

    return fits_file


def masked_name(filename):
    """
    Return foo_masked.fits from foo.fits.
    """
    root, ext = os.path.splitext(filename)
    return root + '_masked' + ext


# ============================================================
# Prepare TP cube
# ============================================================

def prep_tp_cube(TP_fits=None, int_cube=None, int_pb=None):
    """
    Prepare the TP cube for feathering.

    Parameters
    ----------
    TP_fits : str
        Original TP FITS cube.
    int_cube : str
        CASA interferometer image used as regridding template.
    int_pb : str
        CASA interferometer PB image.

    Steps
    -----
    1. Import TP FITS -> CASA image
    2. Regrid to interferometer cube
    3. Reorder axes to match PB
    4. Multiply by PB
    5. Export all CASA products to FITS

    Returns
    -------
    dict
        CASA and FITS filenames for generated products.
    """

    root, ext = os.path.splitext(TP_fits)

    TP_image = root + '.image'
    regrid_image = root + '_regrid.image'
    regrid_trans_image = root + '_regrid_trans.image'
    multiplied_image = root + '_regrid_multiplied.image'

    # --------------------------------------------------------
    # Import fits files
    # --------------------------------------------------------

    importfits(
        fitsimage=TP_fits,
        imagename=TP_image,
        overwrite=True
    )

    # --------------------------------------------------------
    # Regrid TP cube to interferometer cube
    # --------------------------------------------------------

    imregrid(
        imagename=TP_image,
        template=int_cube,
        output=regrid_image,
        overwrite=True
    )

    # --------------------------------------------------------
    # Reorder axes
    #
    # RA, Dec, Freq, Stokes
    # ->
    # RA, Dec, Stokes, Freq
    # --------------------------------------------------------

    imtrans(
        imagename=regrid_image,
        outfile=regrid_trans_image,
        order='0132'
    )

    # --------------------------------------------------------
    # Multiply reordered TP cube by interferometer pb
    # --------------------------------------------------------

    immath(
        imagename=[regrid_trans_image, int_pb],
        expr='IM0*IM1',
        outfile=multiplied_image
    )

    # --------------------------------------------------------
    # Export all casa images to fits
    # --------------------------------------------------------

    TP_image_fits = casa_to_fits(TP_image)
    regrid_fits = casa_to_fits(regrid_image)
    regrid_trans_fits = casa_to_fits(regrid_trans_image)
    multiplied_fits = casa_to_fits(multiplied_image)

    return {
        'tp_image': TP_image,
        'tp_image_fits': TP_image_fits,

        'regrid_image': regrid_image,
        'regrid_fits': regrid_fits,

        'regrid_trans_image': regrid_trans_image,
        'regrid_trans_fits': regrid_trans_fits,

        'multiplied_image': multiplied_image,
        'multiplied_fits': multiplied_fits,
    }


# ============================================================
# Feathering
# ============================================================

def feather_cubes(
    TP_cube_regrid_mult=None,
    int_cube=None,
    feathered_cube=None,
    smoothed_cube=None
):
    """
    Feather the TP cube to the interferometer-only cube and
    smooth the result to the TP resolution.

    All input/output image arguments here are CASA images.

    Returns
    -------
    dict
        CASA and FITS filenames for feathered products.
    """

    # --------------------------------------------------------
    # Feather
    # --------------------------------------------------------

    feather(
        imagename=feathered_cube,
        highres=int_cube,
        lowres=TP_cube_regrid_mult,
        sdfactor=1.0
    )

    # --------------------------------------------------------
    # Get TP beam
    # --------------------------------------------------------

    TP_beam = imhead(TP_cube_regrid_mult)['restoringbeam']

    # --------------------------------------------------------
    # Smooth feathered cube to TP resolution
    # --------------------------------------------------------
    
    imsmooth(
        imagename=feathered_cube,
        major=str(TP_beam['major']['value']) + 'arcsec',
        minor=str(TP_beam['minor']['value']) + 'arcsec',
        pa=str(TP_beam['positionangle']['value']) + 'deg',
        targetres=True,
        outfile=smoothed_cube
    )
    
    # --------------------------------------------------------
    # Export casa images to fits
    # --------------------------------------------------------

    feathered_fits = casa_to_fits(feathered_cube)
    smoothed_fits = casa_to_fits(smoothed_cube)

    return {
        'feathered_image': feathered_cube,
        'feathered_fits': feathered_fits,

        'smoothed_image': smoothed_cube,
        'smoothed_fits': smoothed_fits,
    }


# ============================================================
# Calculate interferometer image rms
# ============================================================

def calc_rms(int_cube):
    """
    Calculate RMS from CASA interferometer cube. Use emission-free channels.
    """

    rms = imstat(int_cube, chans='5~10')['rms'][0]

    return rms


# ============================================================
# Assessment mask
# ============================================================

def assessment_mask(
    TP_cube=None,
    int_cube=None,
    smoothed_feathered_image=None,
    int_pb=None,
    rms_int=None
):
    """
    Create a mask defining where feathering should be assessed.

    All file arguments to this function are FITS files.

    Threshold:
        3 * array_rms * (1 / array_pb)
          * (TP beam area / array beam area)

    Returns
    -------
    masked_tp_image : str
        FITS filename.
    masked_feathered_image : str
        FITS filename.
    """

    # --------------------------------------------------------
    # Read FITS data
    # --------------------------------------------------------

    pb_data = fits.getdata(int_pb)

    tp_data, tp_header = fits.getdata(
        TP_cube,
        header=True
    )

    feathered_data, feathered_header = fits.getdata(
        smoothed_feathered_image,
        header=True
    )

    int_header = fits.getheader(int_cube)

    # --------------------------------------------------------
    # Beam areas
    # --------------------------------------------------------

    bmaj_int = int_header['BMAJ']
    bmin_int = int_header['BMIN']

    bmaj_tp = tp_header['BMAJ']
    bmin_tp = tp_header['BMIN']

    beam_int = bmaj_int * bmin_int
    beam_tp = bmaj_tp * bmin_tp

    # --------------------------------------------------------
    # Threshold
    # --------------------------------------------------------

    thresh = (
        3.0
        * rms_int
        * (1.0 / pb_data)
        * (beam_tp / beam_int)
    )

    # --------------------------------------------------------
    # Assessment mask
    # --------------------------------------------------------

    mask = (
        (tp_data > thresh)
        & (feathered_data > thresh)
    )

    masked_tp_data = np.where(mask, tp_data, np.nan)
    masked_feathered_data = np.where(
        mask,
        feathered_data,
        np.nan
    )

    # --------------------------------------------------------
    # Output filenames
    # --------------------------------------------------------

    masked_tp_image = masked_name(TP_cube)

    masked_feathered_image = masked_name(
        smoothed_feathered_image
    )

    # --------------------------------------------------------
    # Write FITS files
    # --------------------------------------------------------

    fits.writeto(
        masked_tp_image,
        masked_tp_data,
        tp_header,
        overwrite=True
    )

    fits.writeto(
        masked_feathered_image,
        masked_feathered_data,
        feathered_header,
        overwrite=True
    )

    return masked_tp_image, masked_feathered_image


# ============================================================
# Accuracy and fidelity
# ============================================================

def accuracy_fidelity_maps(
    smoothed_masked_feathered_image=None,
    masked_tp_image=None,
    line=None
):
    """
    Make accuracy and fidelity parameter cubes.

    Inputs and outputs are FITS files.

    Accuracy:
        A = (I_comb - I_TP) / |I_TP|

    Fidelity:
        F = |I_TP / (I_comb - I_TP)|
    """

    TP, header = fits.getdata(
        masked_tp_image,
        header=True
    )

    comb_smooth = fits.getdata(
        smoothed_masked_feathered_image
    )

    # --------------------------------------------------------
    # Accuracy parameter
    # --------------------------------------------------------

    with np.errstate(divide='ignore', invalid='ignore'):

        accuracy = (
            (comb_smooth - TP)
            / np.abs(TP)
        )

        fidelity = np.abs(
            TP / (comb_smooth - TP)
        )

    # Replace infinities with NaN
    accuracy[~np.isfinite(accuracy)] = np.nan
    fidelity[~np.isfinite(fidelity)] = np.nan

    accuracy_file = (
        f'accuracy_param_masked_{line}.fits'
    )

    fidelity_file = (
        f'fidelity_param_masked_{line}.fits'
    )

    fits.writeto(
        accuracy_file,
        accuracy,
        header,
        overwrite=True
    )

    fits.writeto(
        fidelity_file,
        fidelity,
        header,
        overwrite=True
    )

    return accuracy_file, fidelity_file


# ============================================================
# Accuracy parameter vs. TP flux
# ============================================================

def a_par_vs_flux(
    masked_tp_image=None,
    accuracy_param_image=None,
    line=None,
    nbins=20
):
    """
    Plot accuracy parameter as a function of TP flux.

    Inputs are FITS files.
    """

    cube_tp = fits.getdata(masked_tp_image)
    cube_accuracy = fits.getdata(accuracy_param_image)

    tp_data = cube_tp.flatten()
    accuracy_data = cube_accuracy.flatten()

    mask = (
        np.isfinite(tp_data)
        & np.isfinite(accuracy_data)
        & (tp_data > 0)
    )

    tp_data = tp_data[mask]
    accuracy_data = accuracy_data[mask]

    # --------------------------------------------------------
    # Logarithmic TP-flux bins
    # --------------------------------------------------------

    bin_edges = np.logspace(
        np.log10(tp_data.min()),
        np.log10(tp_data.max()),
        nbins + 1
    )

    tp_means = []
    accuracy_means = []
    accuracy_stds = []

    for i in range(nbins):

        bin_mask = (
            (tp_data >= bin_edges[i])
            & (tp_data < bin_edges[i + 1])
        )

        if np.any(bin_mask):

            tp_means.append(
                np.mean(tp_data[bin_mask])
            )

            accuracy_means.append(
                np.mean(accuracy_data[bin_mask])
            )

            accuracy_stds.append(
                np.std(accuracy_data[bin_mask])
            )

        else:

            tp_means.append(np.nan)
            accuracy_means.append(np.nan)
            accuracy_stds.append(np.nan)

    tp_means = np.asarray(tp_means)
    accuracy_means = np.asarray(accuracy_means)
    accuracy_stds = np.asarray(accuracy_stds)

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plt.figure(figsize=(6, 5))

    plt.errorbar(
        tp_means,
        accuracy_means,
        yerr=accuracy_stds,
        fmt='o',
        capsize=3,
        label='Binned averages'
    )

    plt.xscale('log')
    plt.yscale('linear')

    plt.xlabel("TP flux")
    plt.ylabel("A-par")

    plt.grid(
        True,
        which="both",
        ls="--",
        alpha=0.5
    )

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        f'{line}_a_par_v_tp.png'
    )

    plt.close()


# ============================================================
# Feathered versus TP flux
# ============================================================

def feathered_vs_tp_flux(
    masked_tp_image=None,
    masked_smoothed_feathered_image=None,
    line=None,
    nbins=20
):
    """
    Plot smoothed feathered flux versus TP flux.

    Inputs are FITS files.
    """

    cube_tp = fits.getdata(masked_tp_image)

    cube_feathered = fits.getdata(
        masked_smoothed_feathered_image
    )

    tp_data = cube_tp.flatten()
    feathered_data = cube_feathered.flatten()

    mask = (
        np.isfinite(tp_data)
        & np.isfinite(feathered_data)
        & (tp_data > 0)
        & (feathered_data > 0)
    )

    tp_data = tp_data[mask]
    feathered_data = feathered_data[mask]

    # --------------------------------------------------------
    # Logarithmic TP bins
    # --------------------------------------------------------

    bin_edges = np.logspace(
        np.log10(tp_data.min()),
        np.log10(tp_data.max()),
        nbins + 1
    )

    tp_means = []
    feathered_means = []
    feathered_stds = []

    for i in range(nbins):

        bin_mask = (
            (tp_data >= bin_edges[i])
            & (tp_data < bin_edges[i + 1])
        )

        if np.any(bin_mask):

            tp_means.append(
                np.mean(tp_data[bin_mask])
            )

            feathered_means.append(
                np.mean(feathered_data[bin_mask])
            )

            feathered_stds.append(
                np.std(feathered_data[bin_mask])
            )

        else:

            tp_means.append(np.nan)
            feathered_means.append(np.nan)
            feathered_stds.append(np.nan)

    tp_means = np.asarray(tp_means)
    feathered_means = np.asarray(feathered_means)
    feathered_stds = np.asarray(feathered_stds)

    # --------------------------------------------------------
    # Plot
    # --------------------------------------------------------

    plt.figure(figsize=(6, 5))

    plt.errorbar(
        tp_means,
        feathered_means,
        yerr=feathered_stds,
        fmt='o',
        capsize=3,
        label='Binned averages'
    )

    x_min = np.nanmin(tp_means)
    x_max = np.nanmax(tp_means)

    plt.plot(
        [x_min, x_max],
        [x_min, x_max],
        'k--'
    )

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel("TP flux")
    plt.ylabel("Feathered image flux")

    plt.legend()
    plt.tight_layout()

    plt.savefig(
        f'{line}_feather_v_tp.png'
    )

    plt.close()


# ============================================================
# Main loop
# ============================================================

lines = ['CO', '13CO', 'C18O']

tp_fits = {
    'CO': glob.glob(
        'member*spw17*.sd.im.fits'
    )[0],

    '13CO': glob.glob(
        'member*spw25*.sd.im.fits'
    )[0],

    'C18O': glob.glob(
        'member*spw27*.sd.im.fits'
    )[0]
}


for line in lines:

    print(f'\nProcessing {line}')

    # ========================================================
    # Original CASA products
    # ========================================================

    int_cube = (
        f'mosaic1_M17_{line}_cube_I_int_'
        f'cyclefactor1_scales05.image'
    )

    int_pb = (
        f'mosaic1_M17_{line}_cube_I_int_'
        f'cyclefactor1_scales05.pb'
    )

    TP_fits = tp_fits[line]

    # ========================================================
    # Export interferometer products to FITS
    # ========================================================

    int_cube_fits = casa_to_fits(int_cube)

    int_pb_fits = casa_to_fits(
        int_pb,
        fits_file=int_pb + '.fits'
    )

    # ========================================================
    # Prepare TP cube
    # ========================================================

    tp_products = prep_tp_cube(
        TP_fits=TP_fits,
        int_cube=int_cube,
        int_pb=int_pb
    )

    TP_cube_regrid_mult = (
        tp_products['multiplied_image']
    )

    TP_cube_regrid_mult_fits = (
        tp_products['multiplied_fits']
    )

    # ========================================================
    # Feather
    # ========================================================

    feathered_cube = (
        f'mosaic1_M17_{line}_feathered.image'
    )

    smoothed_cube = (
        f'mosaic1_M17_{line}_'
        f'feathered_smoothed.image'
    )

    feather_products = feather_cubes(
        TP_cube_regrid_mult=TP_cube_regrid_mult,
        int_cube=int_cube,
        feathered_cube=feathered_cube,
        smoothed_cube=smoothed_cube
    )

    smoothed_cube_fits = (
        feather_products['smoothed_fits']
    )

    # ========================================================
    # Calculate interferometer RMS using CASA cube
    # ========================================================

    rms = calc_rms(int_cube)

    print(f'{line} RMS = {rms}')

    # ========================================================
    # Assessment mask
    # ========================================================

    (
        masked_tp_image,
        masked_smoothed_feathered_image
    ) = assessment_mask(
        TP_cube=TP_cube_regrid_mult_fits,
        int_cube=int_cube_fits,
        smoothed_feathered_image=smoothed_cube_fits,
        int_pb=int_pb_fits,
        rms_int=rms
    )

    # ========================================================
    # Accuracy and fidelity maps
    # ========================================================

    (
        accuracy_param_image,
        fidelity_param_image
    ) = accuracy_fidelity_maps(
        smoothed_masked_feathered_image=(
            masked_smoothed_feathered_image
        ),
        masked_tp_image=masked_tp_image,
        line=line
    )

    # ========================================================
    # Diagnostic plots
    # ========================================================

    a_par_vs_flux(
        masked_tp_image=masked_tp_image,
        accuracy_param_image=accuracy_param_image,
        line=line
    )

    feathered_vs_tp_flux(
        masked_tp_image=masked_tp_image,
        masked_smoothed_feathered_image=(
            masked_smoothed_feathered_image
        ),
        line=line
    )

    print(f'Finished {line}')