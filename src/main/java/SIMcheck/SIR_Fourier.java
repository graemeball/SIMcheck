/* 
 * Copyright (c) 2013, Graeme Ball and Micron Oxford,
 * University of Oxford, Department of Biochemistry.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see http://www.gnu.org/licenses/ .
 */

package SIMcheck;

import ij.*;
import ij.plugin.*;
import ij.plugin.filter.GaussianBlur;
import ij.process.*;
import ij.measure.*;
import ij.gui.*;

import java.awt.Color;
import java.awt.image.IndexColorModel;
import java.io.*;

/** 
 * This plugin plots FFTs of reconstructed stacks with resolution rings.
 * @author Graeme Ball <graemeball@gmail.com>
 * @see ij.plugin.FFT
 * @see ij.process.FHT
 */
public class SIR_Fourier implements PlugIn, EProcessor {

    String name = "Reconstructed Data Fourier Plots";
    ResultSet results = new ResultSet(name);
    double[] resolutions = {0.10, 0.12, 0.15, 0.2, 0.3, 0.6};
    double blurRadius = 6.0d;  // default for 512x512
    static final String[] setMinChoices = {"mode", "mean", "min"};
    int setMinChoice = 0;
    
    @Override
    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        GenericDialog gd = new GenericDialog("SIR_Fourier");
        imp.getWidth();
        gd.addNumericField("Blur radius (512x512)", blurRadius, 1);
        gd.addChoice("Crop minimum to", setMinChoices, "mode");
        gd.showDialog();
        if (gd.wasOKed()) {
            this.blurRadius = gd.getNextNumber();
            this.setMinChoice = gd.getNextChoiceIndex();
	        results = exec(imp);
	        results.report();
        }
    }

    /** 
     * Execute plugin functionality: for SIR and ortho-resliced, perform 
     * 2D FFT on each slice, blur, apply LUT and draw resolution rings.
     * @param imps reconstructed SIR data ImagePlus should be first imp
     * @return ResultSet containing FFT imp, ortho FFT imp, radial profile plot 
     */
    public ResultSet exec(ImagePlus... imps) {
        Calibration cal = imps[0].getCalibration();
        ImagePlus imp2 = imps[0].duplicate();
        new StackConverter(imp2).convertToGray32();  // for OrthoReslicer
        
        /// 1. reslice for orthogonal view
        OrthoReslicer orthoReslicer = new OrthoReslicer();
        ImagePlus impOrtho = imp2.duplicate();
        impOrtho = orthoReslicer.exec(impOrtho, false);
        Calibration calOrtho = impOrtho.getCalibration();

        /// 2. FFT slices in SIR and Z-resliced SIR images
        IJ.showStatus("Fourier transforming slices (lateral view)");
        ImagePlus impF = FFT2D.fftImp(imp2);
        IJ.showStatus("Fourier transforming slices (orthogonal view)");
        ImagePlus impOrthoF = FFT2D.fftImp(impOrtho);
        
        /// 3. optionally blur slices and set display range
        blurRadius *= (double)impF.getWidth() / 512.0d;
        IJ.showStatus("Blurring & rescaling slices (lateral view)");
        impF = displaySettings(impF);  // FIXME, scale info
        IJ.showStatus("Blurring & rescaling slices (orthogonal view)");
        impOrthoF = displaySettings(impOrthoF);
        
        /// 4. annotate Fourier transformed slices with resolution rings
        impOrthoF = resizeAndPad(impOrthoF, cal);
        calOrtho.pixelHeight = calOrtho.pixelWidth;  // after resizeAndPad
        impF = overlayResRings(impF, cal);
        impOrthoF = overlayResRings(impOrthoF, calOrtho);
        I1l.copyStackDims(imps[0], impF);
        I1l.copyStackDims(imps[0], impOrthoF);
        results.addInfo(
                "Fourier plots (XY and XZ)", 
                " to estimate resolution & check for artifacts\n"
                + "  - concentric, evenly-spaced colors indicate high resolution\n"
                + "  - flat Fourier spectrum (broad middle rings) indicate weaker\n"
                + "    high frequency information and poor resolution\n"
                + "  - spots in XY Fourier spectrum indicate periodic XY patterns\n"
                + "  - asymmetric FFT indicates decreased resolution due to: angle to angle intensity\n"
                + "    variations, angle-specific k0 error, or angle-specific z-modulation issues\n");
        String fourierLUTfile = "SIMcheckFourier.lut";
        if (blurRadius > 0) {
            IndexColorModel LUT = loadLut(fourierLUTfile);
            double[] displayRange = {0.0d, 255.0d};  // show all
            I1l.applyLUT(impF, LUT, displayRange);
            I1l.applyLUT(impOrthoF, LUT, displayRange);
        }
        impF.setTitle(I1l.makeTitle(imps[0], "FTL"));
        results.addImp("lateral (XY)", impF);
        impOrthoF.setTitle(I1l.makeTitle(imps[0], "FTO"));
        results.addImp("orthogonal / axial (XZ)", impOrthoF);
        ImagePlus radialProfiles = makeRadialProfiles(impF);
        radialProfiles.setTitle(I1l.makeTitle(imps[0], "FTR"));
        results.addImp("lateral Fourier radial profiles (central Z)", 
                radialProfiles);
        return results;
    }
    
    /** Make radial profile plot for each channel at central Z. */
    ImagePlus makeRadialProfiles(ImagePlus imp) {
        int nc = imp.getNChannels();
        int nz = imp.getNSlices();
        ImagePlus[] profiles = new ImagePlus[nc];
        Radial_Profile radialProfiler = new Radial_Profile();
        for (int c = 1; c <= nc; c++) {
            imp.setPosition(c, nz / 2, 1);
            ImageProcessor ip = imp.getProcessor();
            ImagePlus impC = new ImagePlus("impC" + c, ip);
            impC.setProperty("FHT", "F");  // tell radialProfiler it's Fourier
            I1l.copyCal(imp, impC);
            profiles[c - 1] = radialProfiler.exec(impC);
        }
        ImagePlus impProfiles = I1l.mergeChannels(
                "radial profiles", profiles);
        return impProfiles;
    }
    
    /** Resize impOrthoF for same y pixel size as impF & pad square with 0s */
    ImagePlus resizeAndPad(ImagePlus impOrthoF, Calibration cal) {
        int width = impOrthoF.getHeight();
        int height = impOrthoF.getHeight();
        int depth = impOrthoF.getNSlices();
        double rescaleFactor = cal.pixelHeight / cal.pixelDepth;
        int rescaledHeight = (int)((double)height * rescaleFactor);
        IJ.run(impOrthoF, "Scale...", 
                "x=1.0 y=" + rescaleFactor 
                + " z=1.0 width=" + width 
                + " height=" + rescaledHeight
                + " depth=" + depth + " interpolation=Bilinear"
                + " average process create title=impOrthoResized");
        ImagePlus impOrthoFrsz = IJ.getImage();
        impOrthoFrsz.hide();
        int slices = impOrthoFrsz.getStackSize();
        ImageStack stack = impOrthoFrsz.getStack();
        ImageStack padStack = new ImageStack(width, height, 
        		impOrthoF.getStackSize());
        for (int s = 1; s <= slices; s++) {
            ImageProcessor ip = stack.getProcessor(s);
            int insertStart = width * (((height - rescaledHeight) / 2) - 1);
            int insertEnd = insertStart + width * rescaledHeight;
            ImageProcessor pip = new ByteProcessor(width, height);  // to pad
            byte[] pix = (byte[])((ByteProcessor)ip).getPixels();
            byte[] padpix = new byte[width * height];
            for (int i = insertStart; i < insertEnd; i++) {
                padpix[i] = pix[i - insertStart];
            }
            pip.setPixels(padpix);
            padStack.setProcessor(pip, s);
        }
        impOrthoFrsz.setStack(padStack);
        I1l.copyCal(impOrthoF, impOrthoFrsz);
        impOrthoFrsz.setProperty("FHT", impOrthoF.getProperty("FHT"));
        return impOrthoFrsz;
    }

    /** Optional gaussian blur, and select lower end of intensity range. */
    ImagePlus displaySettings(ImagePlus imp) {
        int ns = imp.getStackSize();
        GaussianBlur gblur = new GaussianBlur();
        gblur.showProgress(false);
        for (int s = 1; s <= ns; s++) {
            imp.setSlice(s);
            if (blurRadius > 0) {
                ImageProcessor ip = imp.getProcessor();
                gblur.blurGaussian(ip, blurRadius, blurRadius, 0.002);
            }
            ImageStatistics stats = imp.getProcessor().getStatistics();
            double min;
            switch (setMinChoice) {
                case 0:  min = (double)stats.mode;
                         break;
                case 1:  min = (double)stats.mean;
                         break;
                case 2:  min = (double)stats.min;
                         break;
                default: min = (double)stats.min;
                         break;
            }
            double max = (double)stats.max;
            ByteProcessor bp = (ByteProcessor)imp.getProcessor();
            ImageProcessor ip = 
                    (ImageProcessor)setBPminMax(bp, (int)min, (int)max, 255);
            imp.setProcessor(ip);
            IJ.showProgress(s, ns);
        }
        return imp;
    }
    
    /** Set ByteProcessor range min to inMax to output range 0 to outMax. */
    ByteProcessor setBPminMax(ByteProcessor bp, 
            int min, int inMax, int outMax) {
        if (min < 0 || inMax > 255 || outMax > 255) {
            throw new IllegalArgumentException("invalid min or max for 8-bit");
        }
        byte[] bpix = (byte[])bp.getPixels();
        int range = inMax - min; 
        for (int i = 0; i < bpix.length; i++) {
            int scaledPix = (int)bpix[i] & 0xff;
            if (scaledPix > inMax) {
                scaledPix = inMax;
            }
            scaledPix -= min;
            if (scaledPix < 0) {
                scaledPix = 0;
            }
            scaledPix = (int)(outMax * ((double)scaledPix / range));
            bpix[i] = (byte)scaledPix;
        }
        bp.setPixels(bpix);
        return bp;
    }
  
    /** 
     * Overlay resolution rings on each slice of a Fourier ImagePlus. 
     *  NB. raw (non-FFT) imp is required for original calibrations.
     */
    ImagePlus overlayResRings(ImagePlus Fimp, Calibration cal) {
        String unit = cal.getUnit();
        if (!(unit.startsWith(""+IJ.micronSymbol) || unit.startsWith("u") 
                || unit.startsWith("micro"))) {
            IJ.log("  ! warning - non-micron calibration (" + unit
                    + ") - cannot plot resolutions");
        } else {
            int width = Fimp.getWidth();
            int height = Fimp.getHeight();
            double pixWidth = cal.pixelWidth;
            double pixHeight = cal.pixelHeight;
            int fontSize = Fimp.getProcessor().getFont().getSize();
            Overlay resOverlay = new Overlay();
            for (int ring = 0; ring < resolutions.length; ring++) {
                // res = pixel size * width/radius (i.e. 2cycles)
                double dresX = ((double) width * pixWidth 
                		/ resolutions[ring]);
                int resX = (int) dresX;
                double dresY = ((double) height * pixHeight 
                		/ resolutions[ring]);
                int resY = (int) dresY;
                // NB. referring to *display* X & Y - could be e.g. true Z
                int ovalW = resX * 2;
                int ovalH = resY * 2;
                int topLeftX = width / 2 + 1 - ovalW / 2;
                int topLeftY = height / 2 + 1 - ovalH / 2;
                OvalRoi currentOval = new OvalRoi(topLeftX, topLeftY, 
                        ovalW, ovalH);
                currentOval.setStrokeColor(Color.WHITE);
                currentOval.setStrokeWidth(1.0);
                resOverlay.add(currentOval);
                TextRoi currentRes = new TextRoi(
                        (int) (width / 2 - fontSize),
                        (int) (height / 2 - 2 * (-0.5 + (ring + 1) % 2)
                                * (ovalH / 2) - fontSize),
                        Double.toString(resolutions[ring]));
                currentRes.setStrokeColor(Color.WHITE);
                resOverlay.add(currentRes);
            }
            Fimp.setOverlay(resOverlay);
        }
        return Fimp;
    }

    /** Load a LUT from a file. (NB. getClass is non-static) */
    IndexColorModel loadLut(String LUTfile) {
        InputStream is = getClass().getResourceAsStream(LUTfile);
        IndexColorModel cm = null;
        if (is != null) {
            try {
                cm = LutLoader.open(is);
            } catch (IOException e) {
                IJ.log("  ! error loading LUT");
                IJ.error("" + e);
            }
        }
        if (cm == null) {
            IJ.log("  ! error loading LUT");
        }
        return cm;
    }
    
    /** main method for testing */
    public static void main(String[] args) {
        System.out.println("Testing SIR_Fourier.java");
        new ImageJ();
        ImagePlus impTest = IJ.openImage(
        		"/Users/graemeb/Workspace/SIMcheck/test_images/Test/"
        		+ "V3_DAPI_good_3Dsmall.tif");
        impTest.show();
        OrthoReslicer orthoReslicer = new OrthoReslicer();
        ImagePlus impOrtho = orthoReslicer.exec(impTest, false);
        impOrtho.show();
        IJ.log("impOrtho height/width after reslice: " + impOrtho.getHeight() 
        		+ "/" + impOrtho.getWidth());
        ImagePlus impOrthoF = FFT2D.fftImp(impOrtho.duplicate());
        impOrthoF.show();
        int width = impOrtho.getWidth();
        int height = impOrtho.getHeight();
        ImageStack windowedStack = new ImageStack(width, height);
        int paddedSize = FFT2D.calcPadSize(impOrtho);
        if (paddedSize != width || paddedSize != height) {
            width = paddedSize;
            height = paddedSize;
        }
        ImageStack paddedStack = new ImageStack(width, height);
        for (int s = 1; s <= impOrtho.getStackSize(); s++) {
            ImageProcessor ip = impOrtho.getStack().getProcessor(s);
            ImageProcessor wp = FFT2D.gaussWindow(ip.duplicate(), 0.125d);
            windowedStack.addSlice(wp);
            ImageProcessor pp = FFT2D.pad(wp.duplicate(), paddedSize); 
            paddedStack.addSlice(pp);
        }
        ImagePlus impWindowed = new ImagePlus();
        impWindowed.setStack(windowedStack);
        impWindowed.setTitle("windowed");
        impWindowed.show();
        ImagePlus impPadded = new ImagePlus();
        impPadded.setStack(paddedStack);
        impPadded.setTitle("padded");
        impPadded.show();
    }
}