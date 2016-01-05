# Output files #
  * image\_ref.raw - Reference image of emitted decays that were detected in the coincidences of the ROOT file (Not reconstructed).
  * image\_ref\_sens\_corr.raw - Same as before, but corrected by sensitivity
  * image.raw - Reconstructed image (without sensitivity correction)
  * image\_norm.raw - Reconstructed image with sensitivity correction
  * image\_norm\_med.raw - Reconstructed image with sensitivity correction with median filter applied at the end.

# Format #
  * All images are 75x75x61 voxels, 32bits, little endian.
  * All images are float, except for image.raw (integer)

# Tool for MATLAB #

> The tool fbi3d\_raw2mat.m allows converting these images into a .mat file for MATLAB