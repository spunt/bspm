INSTALLATION

This simple but useful utility toolbox consists of three files:
- README.txt:
  this file;
- tbx_cfg_logtransform.m:
  a configuration function;
- tbx_run_logtransform.m:
  a job execution function.
Put them in a subfolder named "Logtransform" of the "toolbox" folder in the SPM8/12 distribution. After starting spm, an additional menu item called "Log-transform" will appear in the SPM > Tools menu of the Batch Editor window. You may include this transformation as a module in your batch processing pipeline.


FUNCTION

The function of this utility is to apply a logarithmic transformation to fMRI image volumes according to Y = 100*ln(X/X0), where X is a voxel's original image intensity, X0 is a scaling factor, and Y is the voxel's output intensity. This transformation converts changes relative to a baseline into absolute changes. In other words, whereas the input images are expressed in arbitrary scanner units, the output images are in units of percentage signal change.
To see this, suppose a voxel has a baseline signal B on top of which a BOLD response signal R is imposed. Then
Y = 100 * ln( X / X0 )
  = 100 * ln( (B+R) / X0 )
  = 100 * ln( B/X0 * (1+R/B) )
  = 100 * ln(B/X0) + 100 * ln(1+R/B)
  ~ 100 * ln(B/X0) + 100 * R/B
The latter approximation can be made because BOLD signal fluctuations R are much smaller than the baseline signal B itself.
The resulting signal Y can be seen to consist of a transformed baseline 100*ln(B/X0), which will be accounted for in the GLM, plus a transformed response 100*R/B, which is precisely the response R expressed as a percentage relative to the baseline B.
Therefore, the response-related regression coefficients beta that are estimated in the GLM can directly be interpreted as percentage signal change, without further need for an explicit division by the baseline signal. These values can directly be compared across subjects and imported to second-level analyses. Thus, the conversion to percentage signal change can be taken care of during preprocessing already, simplifying the analysis pipeline.


USAGE

A Log-transform module will contain the following items:
- Data:
  a structure containing an arbitrary number of Session subitems, each of which contain lists of image volumes (the specification of input images is analogous to that in a Spatial > Realign module);
- Scaling:
  The scaling factor X0 can either be determined automatically by the algorithm or be specified by the user, either as a constant value or as an image with voxel-wise scaling factors;
- Negative values:
  After the transformation, negative values may occur that can either be retained or clipped to zero;
- Data type
  The data type of the output images can either be kept the same as the corresponding input images, or be set to one of the predefined SPM data types.
Both manual file selection and dependencies between modules are supported.


NOTES

- In principle, the value of the scaling factor X0 is irrelevant. Because Y = 100*ln(X/X0) can also be written as Y = 100*ln(X)-100*ln(X0), the scaling factor can be seen to act as an additive offset. Such terms are normally incorporated into the baseline of the GLM, which is not of interest given that the responses are already expressed in units of percentage signal change after the log-transform.
- Despite the aforementioned, a proper choice of the scaling factor X0 should be made for practical reasons. First, this may help guarantee that small responses are not lost due to discretisation and round-off. Second, after the log-transformation, image intensities may attain negative values, that are not always well handled by third-party software. As a rule of thumb, the scaling factor X0 should be chosen to lie somewhat below the input intensity of the grey matter segment. This will result in transformed values that are positive but not too far from zero. The default Automatic scaling option should normally be able to achieve this.
- When using Automatic scaling, it is important that all images in a single session have similar gain and contrast, because a suitable scaling factor X0 will be estimated from the first image and subsequently applied to all images in a session. At the same time, images that belong to the same imaging run should not be split up over sessions, because otherwise different scaling factors will be applied to them, rendering their direct comparison meaningless. When scaling factors are specified by the user, the same transformation is applied to all sessions, so it makes no difference whether images are spread over multiple session or combined in a single one.
- Other sources of signal fluctuations, like scanner drift or residual motion, should be dealt with in the GLM as usual, for instance by including them as covariates or by filtering the voxel signals. The use of the transformation should not be combined with other methods to obtain percentage signal change measures, like Global normalisation to a mean of 100.
- The data can be written in a number of different SPM file formats, including the same format as the input images. The 8-bit formats UINT8 and INT8 will lead to large discretisation errors and may not be able to represent response magnitudes sufficiently accurately. The use of the default SAME option will usually be appropriate.
- The output of negative transformed voxel intensities can be suppressed by clipping them to zero using the Negative values option. Using a proper scaling factor X0, the output images may thus retain a natural tissue contrast. Use clipping if negative values are improperly dealt with or when uninformative voxels need to be excluded. Alternatively, choose UINT8, UINT16, or UINT32 as the output format, which do not support negative values (these will also allow the sign bit to be used to represent the signal, reducing the discretisation error). Note that the INT8, UINT16, and UINT32 formats are supported by SPM, but third-party software may wrongly interpret the sign bit and display these as UINT8, INT16, and INT32, respectively.
- When using multi-image volumes in combination with Scaling by image, it is a bad idea to scale by the first image of the Session. The first output image would then be a null image, resulting in undefined intensity scaling, prohibiting subsequent images to be properly appended to the same file. Warnings may also occur for multi-volume volumes when the file already exists; the volumes should still be overwritten correctly.


CONTACT

For comments or feedback, feel free to contact me by e-mail: Dave Langers <dave@langers.nl>.
