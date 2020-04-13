# NeXLMatrixCorrection

Implements the XPP matrix correction and Reed fluorescence correction algorithms for bulk and coated samples.

Primarily these algorithms are designed to take a 'Vector{KRatio}' and return a 'Material'.  Since they are
intended for both WDS and EDS, the k-ratio can represent one or more characteristic X-ray lines from a single
element.  K-ratios compare a measured intensity with the intensity from a reference (standard) material. Typically,
these two materials are measured at the same beam energy but multiple beam energy measurements are also supported.
