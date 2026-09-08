from __future__ import annotations
import numpy
from astropy import units


class PhotSpec:
    """
    Represent the spectrum of the objects as a set of colors.

    Historically the colours were in a list with the following order:

    Array of colors (10*R8)
                colors(1) : g-x
                colors(2) : r-x
                colors(3) : i-x
                colors(4) : z-x
                colors(5) : u-x
                colors(6) : V-x
                colors(7) : B-x
                colors(8) : R-x
                colors(9) : I-x
                colors(10) : w-x

    """
    OLD_BAND_ORDER = ['g', 'r', 'i', 'z', 'u', 'V', 'B', 'R', 'I', 'w']

    COLORS = dict([(
        'default', dict([('g-g', 0.0 * units.mag),
                         ('r-g', -0.7 * units.mag),
                         ('i-g', -1.2 * units.mag),
                         ('z-g', -1.7 * units.mag),
                         ('u-g', +0.3 * units.mag),
                         ('V-g', +0.5 * units.mag),
                         ('B-g', +0.1 * units.mag),
                         ('R-g', -0.8 * units.mag),
                         ('I-g', -1.2 * units.mag),
                         ('w-g', -0.8 * units.mag),
                         ('J-g', -1.5 * units.mag),
                         ('H-g', -1.6 * units.mag)])),
        ('implanted', dict([('g-g', 0.0 * units.mag),
                            ('r-g', -0.7 * units.mag),
                            ('i-g', -1.2 * units.mag),
                            ('z-g', -1.7 * units.mag),
                            ('u-g', +0.3 * units.mag),
                            ('V-g', +0.5 * units.mag),
                            ('B-g', +0.1  * units.mag),
                            ('R-g', -0.8 * units.mag),
                            ('I-g', -1.2 * units.mag),
                            ('w-g', -0.8 * units.mag),
                            ('J-g', -1.5 * units.mag),
                            ('H-g', -1.6 * units.mag)])),
        ('cold', dict([('g-g', 0.0 * units.mag),
                       ('r-g', -0.7 * units.mag),
                       ('i-g', -1.2 * units.mag),
                       ('z-g', -1.7 * units.mag),
                       ('u-g', +0.3 * units.mag),
                       ('V-g', +0.5 * units.mag),
                       ('B-g', +0.1 * units.mag),
                       ('R-g', -0.8 * units.mag),
                       ('I-g', -1.2 * units.mag),
                       ('w-g', -0.8 * units.mag),
                       ('J-g', -1.5 * units.mag),
                       ('H-g', -1.6 * units.mag)]))])

    def __init__(self, colors=None):
        if colors is None:
            colors = PhotSpec.COLORS
        self.colors = colors

    @property
    def spectral_groups(self) -> list:
        return list(self.colors)

    def orbital_to_spectral_group(self, orbital_group) -> str:
        """
        Find the color component that is the best match for the given orbital component.
        """
        spectral_group = 'default'
        list_of_matching_spectral_groups = numpy.arange(len(self.spectral_groups))[[x in orbital_group
                                                                                    for x in self.spectral_groups]]
        if len(list_of_matching_spectral_groups) > 0:
            spectral_group = self.spectral_groups[list_of_matching_spectral_groups[0]]
        return spectral_group

    def __repr__(self):
        return str(self.colors)

    def __call__(self, orbital_group, model_band) -> dict:
        """
        Find the colour component that is the best match for the given component.
        """
        return self.transform_spectral_group_to_model_band(self.orbital_to_spectral_group(orbital_group), model_band)

    def transform_spectral_group_to_model_band(self, spectral_group, model_band):
        """
        Set the base bandpass for specphot dictionary of colors
        """
        band_ratios = list(self.colors[spectral_group])
        specphot = {}
        for band_ratio in band_ratios:
            this_band, base_band = band_ratio.split('-')
            this_color = self.colors[spectral_group][band_ratio]
            offset_color = self.colors[spectral_group][f"{model_band}-{base_band}"]
            specphot[f"{this_band}-{model_band}"] = this_color - offset_color
        return specphot

    def colors_list(self, orbital_group: str, model_band: str) -> numpy.array:
        """
        Return the list of colors for the given orbital group where the index of the color value in the list
        is ascii code of the first character of the band-ratio of the color.

        This is used to create the list of colors that is passed into the Fortran component of SSim.
        Layout matches Fortran ``filter_to_index``: Python index ``ord(L)-ord('A')`` is Fortran
        index ``IACHAR(L)-IACHAR('A')+1`` when the list is passed through f2py.
        """
        spec_phot_list = numpy.zeros(ord('z')-ord('A')+1) * units.mag
        colors = self.transform_spectral_group_to_model_band(self.orbital_to_spectral_group(orbital_group), model_band)
        for band_ratio in colors:
            bandpass = band_ratio.split('-')[0]
            if len(bandpass) != 1 or not ('A' <= bandpass <= 'z'):
                raise ValueError(f"Bandpass {bandpass!r} is out of range A-z for the color array")
            spec_phot_list[ord(bandpass) - ord('A')] = colors[band_ratio]
        return list(spec_phot_list.to('mag').value)

    @classmethod
    def from_list(cls, colors_list: numpy.array, model_band: str):
        """
        Rebuild a PhotSpec from a sparse ASCII-indexed color list (as produced by ``colors_list``).

        List index ``i`` corresponds to bandpass ``chr(i + ord('A'))``; non-zero entries become
        ``{band}-{model_band}`` color terms.
        """
        spec_phot = {}
        for idx, color in enumerate(colors_list):
            if color != 0:
                band_ratio = f"{chr(idx + ord('A'))}-{model_band}"
                spec_phot[band_ratio] = color * units.mag if not hasattr(color, 'unit') else color
        return cls(colors={'default': spec_phot})

    @classmethod
    def from_old_style_list(cls, colors_list):
        """
        In previous model files the list of colors was stored as a list of floats and the order of the list
        determined which band the list element referred to.  The default COLORS object is in this historic ordering.

        We pick the smallest absolute color as the base band, yes it might not be zero,
        but what other choice can we use?
        """
        base_band_idx = numpy.argmin(numpy.fabs(colors_list))
        spec_phot = {}
        for idx, color in enumerate(colors_list):
            spec_phot[f"{cls.OLD_BAND_ORDER[idx]}-{cls.OLD_BAND_ORDER[base_band_idx]}"] = color
        return cls(colors={'default': spec_phot})
