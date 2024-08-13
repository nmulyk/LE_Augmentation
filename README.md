The LE_Augmentation repository consists of 3 codes: `manufacture_LE.py`, `augment_real_LE.py`, `mask_bright_stars.py`. These codes were originally used to increase the size of and improve our light echo training set for ALED_TF2, found at https://github.com/AmanjotBhullar/ALED_TF2. For a complete description of this work, please see **insert link to thesis**.

### manufacture_LE.py
`manufacture_LE.py` manufactures Dragonfly Telephoto Array (DTA) light echoes from Canada-France-Hawaii Telescope (CFHT) detections by isolating and overlaying the CFHT light echoes onto DTA images, using the following steps:

1. Isolate the CFHT light echo pixels by using the light echo mask to set non-light echo pixels to zero.
2. Pad the edges of the CFHT image so it is the same shape as the DTA image.
3. Reproject the isolated CFHT light echo to rescale it to match the DTA's scale.
4. Optional: Adjust the brightness of the light echo. This option is used to test whether a dimmer light echo can be recognized by ALED.
5. Optional: Randomly change the light echo's orientation and location. This step helps increase the variation in the manufactured images but, due to the pose invariance property of capsule networks, is not expected to linearly improve training.
6. Change the resolution of the light echo with a Gaussian filter of $\sigma = 1$ to match the DTA image.
7. Overlay the light echo onto the DTA image. As the DTA image only contains the sky and artifacts from differencing bright stars, it serves as the background for the new image.
8. Create a mask for the new light echo by setting all light echo-containing pixels to one and all non-light echo pixels to zero.

 ### augment_real_LE.py
 `augment_real_LE.py` augments the existing real DTA light echo by overlaying the isolated real light echo onto different DTA images while varying the position, orientation, brightness, and width of the gap between linear segments of the light echo, using the following steps:

 1. Create a detailed mask for the real light echo using GIMP (https://www.gimp.org/).
 2. Isolate the real light echo pixels by using the light echo mask to set non-light echo pixels to zero.
 3. Pad the edges of the real light echo image so it is the same shape as the DTA image.
 4. Reproject the isolated real light echo onto the new DTA image.
 5. Optional: Adjust the brightness of the light echo.
 6. Optional: Adjust the width of the gap between line segments of the light echo. This is done by rotating the light echo so that the two light echo line segments are vertical, isolating the pixels associated with each line segment, and shifting the pixels horizontally. By doing this, it mimics a SN with long emission (wide) or short emission (narrow).
 7. Optional: Randomly change the light echo's orientation and location.
 8. Overlay the light echo onto the new DTA image.
 9. Create a mask for the new light echo by setting all light echo-containing pixels to one and all non-light echo pixels to zero.

### mask_bright_stars.py
`mask_bright_stars.py` masks out bright star difference artifacts, which were misidentified as light echoes by ALED, using the following steps:

1. Locate the region encompassed by the image using the image header.
2. Find the same region in the Fourth U.S. Naval Observatory CCD Astrograph Catalog (https://irsa.ipac.caltech.edu/data/UCAC4/ucac4.html) and list the locations of all stars in this region.
3. Create a circular aperture around each star with the radius found with a piece-wise function of the stars' magnitude.
4. Mask out the stars by setting all pixels within the circular apertures to zero.
