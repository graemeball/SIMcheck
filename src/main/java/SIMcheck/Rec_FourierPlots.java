/* 
 * Copyright (c) 2015, Graeme Ball and Micron Oxford,
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
import ij.process.*;
import ij.measure.*;
import ij.gui.*;

import java.awt.Color;
import java.awt.image.IndexColorModel;

/** 
 * This plugin plots FFTs of reconstructed stacks with resolution rings.
 * @author Graeme Ball <graemeball@gmail.com>
 * @see ij.plugin.FFT
 * @see ij.process.FHT
 */
public class Rec_FourierPlots implements PlugIn, Executable {

    public static final String name = "Reconstructed Fourier Plots";
    public static final String TLA = "FTL";  // Fourier Transform Lateral
    public static final String TLA2 = "FTR";  // FT Radial profile
    public static final String TLA3 = "FTO";  // FT Orthogonal (XZ)
    private ResultSet results = new ResultSet(name, TLA);
    private static final IndexColorModel fourierLUT = 
            I1l.loadLut("SIMcheck/SIMcheckFourier.lut");
    
    // parameter fields
    public double[] resolutions = {0.10, 0.12, 0.15, 0.2, 0.3, 0.6};
    public double blurRadius = 6.0d;  // default for 512x512
    public double winFraction = 0.06d;  // window function size, 0-1
    // TODO: refactor & remove default winFraction from here
    
    // options
    public boolean fft3d = true;          // use 3D FFT? (requires ParallelFFTJ)
    public boolean autoCutoff = true;     // no noise cut-off?
    public boolean manualCutoff = false;  // manual noise cut-off?
    public boolean applyWinFunc = false;  // apply window function?
    public boolean blurAndLUT = false;    // blur & apply false color LUT?
    public boolean showAxial = false;     // show axial FFT?
    
    private double[] channelMinima = null;
    
    @Override
    public void run(String arg) {
        ImagePlus imp = IJ.getImage();
        imp.getWidth();
        try {
            Class.forName("edu.emory.mathcs.parallelfftj.FloatTransformer");
        } catch(ClassNotFoundException e) {
            this.fft3d = false;
        }
        boolean doCheck = false;
        GenericDialog gd = new GenericDialog(name);
        gd.addCheckbox("(1)_Cut-off:_auto (stack mode)", autoCutoff);
        gd.addCheckbox("     Cut-off:_manual (default=0)", manualCutoff);
        gd.addCheckbox("(2)_Window_function*", applyWinFunc);
        gd.addCheckbox("(3)_Blur_&_false-color_LUT", blurAndLUT);
        gd.addCheckbox("(4)_Show_axial_FFT", showAxial);
        boolean showedFFT3Doption = false;
        if (fft3d) {
            gd.addCheckbox("(5)_3D_FFT (ParallelFFTJ)", fft3d);
            showedFFT3Doption = true;
        }
        gd.addMessage("* suppresses edge artifacts");
        gd.showDialog();
        if (gd.wasOKed()) {
            doCheck = true;
            // TODO: notCutoff, manualCutoff and autoScale radioButton group
            this.autoCutoff = gd.getNextBoolean();
            this.manualCutoff = gd.getNextBoolean();
            this.applyWinFunc = gd.getNextBoolean();
            this.blurAndLUT = gd.getNextBoolean();
            this.showAxial = gd.getNextBoolean();
            if (showedFFT3Doption) {
                this.fft3d = gd.getNextBoolean();
            }
            if (manualCutoff) {
                // skip if noCutoff
                this.channelMinima = new double[imp.getNChannels()];
                SIMcheck_.specifyBackgrounds(
                        channelMinima, "Set noise cut-off:");
            }
        } else {
            doCheck = false;
            
        }
        if (doCheck) {
            results = exec(imp);
            results.report();
        }
    }
    

    /** 
     * Execute plugin functionality: for reconstructed and ortho-resliced,
     * 2D FFT each slice and draw resolution rings (optional blur+LUT).
     * @param imps reconstructed data ImagePlus should be first imp
     * @return ResultSet containing FFT imp, ortho FFT imp, radial profile plot 
     */
    public ResultSet exec(ImagePlus... imps) {
        try {
            Class.forName("edu.emory.mathcs.parallelfftj.FloatTransformer");
            // use parallel fftj for 3d transform if available
        } catch(ClassNotFoundException e) {
            this.fft3d = false;  // can't find ImgLib2, so no 3D FFT
        }
        // TODO: check we have micron calibrations before continuing..
        Calibration cal = imps[0].getCalibration();
        ImagePlus imp2 = null;
        if (manualCutoff) {
            imp2 = Util_RescaleTo16bit.exec(imps[0].duplicate(), channelMinima);
        } else if (autoCutoff) {
            imp2 = Util_RescaleTo16bit.exec(imps[0].duplicate());
        } else {
            imp2 = imps[0].duplicate();
            IJ.run("Conversions...", "scale");
            IJ.run(imp2, "16-bit", "");
            IJ.run("Conversions...", " ");
        }
        IJ.showStatus("Fourier transforming z-sections (lateral view)");
        if (!applyWinFunc) {
            winFraction = 0.0d;
        }
        ImagePlus impF = null;
        if (fft3d) {
            IJ.run(imp2, "32-bit", "");
            impF = FFT3D.fftImp(imp2, winFraction);
            impF = gaussBlur(impF);
            IJ.setMinAndMax(impF, 20, 40);
            setLUT(impF, 20.0d, 40.0d);
        } else {
            impF = FFT2D.fftImp(imp2, true, winFraction, 0.2d);
            impF = gaussBlur(impF);
            if (imps[0].isComposite()) {
                impF = new CompositeImage(impF);
            }
            IJ.setMinAndMax(impF, 2, 40);
            setLUT(impF, 2.0d, 40.0d);
        }
        impF = overlayResRings(impF, cal);
        I1l.copyStackDims(imps[0], impF);
        impF.setPosition(1, impF.getNSlices() / 2, 1);
        impF.setTitle(I1l.makeTitle(imps[0], TLA));
        results.addImp("Fourier Transform Lateral (XY; resolution rings" +
                " in microns)", impF);
        // radial & axial profiles only if we have calibration info
        if (imp2.getCalibration().getUnit().equals("pixel")) {
            IJ.log("Calibration info required for Radial & Axial FFT plots!");
        } else {
            ImagePlus radialProfiles = makeRadialProfiles(impF);
            radialProfiles.setTitle(I1l.makeTitle(imps[0], TLA2));
            results.addImp("Fourier Transform Radial profile "
                    + "(lateral, central Z)", radialProfiles);
            if (showAxial) {
                Calibration calOrtho = null;
                ImagePlus impOrtho = null;
                ImagePlus impOrthoF = null;
                if (fft3d) {
                    new StackConverter(impF).convertToGray32();  // for OrthoReslicer
                    OrthoReslicer orthoReslicer = new OrthoReslicer();
                    ImagePlus impF2 = impF.duplicate();
                    //impOrthoF = I1l.takeCentralZ(impOrthoF);
                    impOrthoF = orthoReslicer.exec(impF2, true);
                    impOrthoF = resizeAndPad3d(impOrthoF, cal);
                    impOrthoF = gaussBlur(impOrthoF);
                    IJ.setMinAndMax(impOrthoF, 20, 40);
                    setLUT(impOrthoF, 20.0d, 40.0d);
                    calOrtho = impOrthoF.getCalibration();
                } else {
                    /// for orthogonal (axial) view, reslice first
                    new StackConverter(imp2).convertToGray32();  // for OrthoReslicer
                    OrthoReslicer orthoReslicer = new OrthoReslicer();
                    impOrtho = imp2.duplicate();
                    impOrtho = orthoReslicer.exec(impOrtho, false);
                    impOrtho = I1l.takeCentralZ(impOrtho);
                    calOrtho = impOrtho.getCalibration();
                    IJ.showStatus("FFT z-sections (orthogonal view)");
                    impOrthoF = FFT2D.fftImp(impOrtho, true, winFraction, 0.2d);
                    impOrthoF = resizeAndPad(impOrthoF, cal);
                    impOrthoF = gaussBlur(impOrthoF);
                    IJ.setMinAndMax(impOrthoF, 2, 40);
                    setLUT(impOrthoF, 2.0d, 40.0d);
                }
                calOrtho.pixelHeight = calOrtho.pixelWidth;  // after resizeAndPad
                impOrthoF = overlayResRings(impOrthoF, calOrtho);
                I1l.copyStackDims(imps[0], impOrthoF);
                impOrthoF.setPosition(1, impF.getNSlices() / 2, 1);
                impOrthoF.setTitle(I1l.makeTitle(imps[0], TLA3));
                results.addImp("Fourier Transform Orthogonal (XZ)", impOrthoF);
            }
        }
        impF.setPosition(1, impF.getNSlices() / 2, 1);
        if (fft3d) {
            results.addInfo("About",
                    "by default the reconstructed data are (1) cropped to mode;"
                            + " (2) 3D Fourier transformed (ParallelFFTJ)"
                            + " and log-scaled power spectrum displayed.");
        } else {
            results.addInfo("About",
                    "by default the reconstructed data are (1) cropped to mode;"
                            + " (2) Fourier transformed and log-scaled.");
        }
        results.addInfo("How to interpret", 
            "Fourier plots highlight potential artifacts and indicate"
            + " effective resolution:"
            + "  - Spots in Fourier spectrum indicate periodic patterns."
            + "  - Flattened Fourier spectrum (plateau in radial profile)"
            + " indicates lack of high frequency signal and poor resolution."
            + "  - Asymmetric FFT indicates angle-specific decrease in"
            + " resolution due to: angle-to-angle intensity variations,"
            + " angle-specific illumination pattern ('k0') fit error, or"
            + " angle-specific z-modulation issues.  -- ");
        return results;
    }
    
    /** Gaussian-blur result, or not, based on blurAndLUT option field. */
    private ImagePlus gaussBlur(ImagePlus imp) {
        if (blurAndLUT) {
            imp = I1l.gaussBlur2D(imp, blurRadius);
        }
        return imp;
    }
    

    /** Use color LUT, or not; and update display range. */
    private void setLUT(ImagePlus imp, double min, double max) {
        if (blurAndLUT) {
            double[] displayRange = {min, max};
            I1l.applyLUT(imp, fourierLUT, displayRange);
        } else {
            if (imp.isComposite()) {
                CompositeImage ci = (CompositeImage)imp;
                ci.setMode(IJ.GRAYSCALE);
            }
        }
    }

    /** Resize imp for same y pixel size as in cal & pad square with 0s. */
    private ImagePlus resizeAndPad(ImagePlus imp, Calibration cal) {
        int width = imp.getHeight();
        int height = imp.getHeight();
        int depth = imp.getNSlices();
        double rescaleFactor = cal.pixelHeight / cal.pixelDepth;
        int rescaledHeight = (int)((double)height * rescaleFactor);
        IJ.run(imp, "Scale...", 
                "x=1.0 y=" + rescaleFactor 
                + " z=1.0 width=" + width 
                + " height=" + rescaledHeight
                + " depth=" + depth + " interpolation=Bilinear"
                + " average process create title=impOrthoResized");
        ImagePlus imp2 = IJ.getImage();  // TODO: refactor
        imp2.hide();
        int slices = imp2.getStackSize();
        ImageStack stack = imp2.getStack();
        ImageStack padStack = new ImageStack(width, height, 
                imp.getStackSize());
        for (int s = 1; s <= slices; s++) {
            ImageProcessor ip = stack.getProcessor(s);
            int insertStart = width * (((height - rescaledHeight) / 2) + 1);
            int insertEnd = insertStart + width * rescaledHeight;
            if (ip instanceof ByteProcessor) {
                ImageProcessor pip = new ByteProcessor(width, height); // to pad
                byte[] pix = (byte[])((ByteProcessor)ip).getPixels();
                byte[] padpix = new byte[width * height];
                for (int i = insertStart; i < insertEnd; i++) {
                    padpix[i] = pix[i - insertStart];
                }
                pip.setPixels(padpix);
                padStack.setProcessor(pip, s);
            } else {
                ImageProcessor pip = new FloatProcessor(width, height); // to pad
                float[] pix = (float[])((FloatProcessor)ip).getPixels();
                float[] padpix = new float[width * height];
                for (int i = insertStart; i < insertEnd; i++) {
                    padpix[i] = pix[i - insertStart];
                }
                pip.setPixels(padpix);
                padStack.setProcessor(pip, s);
            }
        }
        imp2.setStack(padStack);
        I1l.copyCal(imp, imp2);
        imp2.setProperty("FHT", imp.getProperty("FHT"));
        return imp2;
    }
    
    /** Resize and pad square resliced 3D FFT result. */
    private ImagePlus resizeAndPad3d(ImagePlus imp, Calibration cal) {
        int width = imp.getWidth();
        int height = imp.getHeight();
        int depth = imp.getNSlices();
        double rescaleFactor = cal.pixelHeight / cal.pixelDepth;
        int rescaledHeight = (int)((double)width * rescaleFactor);
        double rescaleY = (double)height / rescaledHeight;
        IJ.run(imp, "Scale...", 
                "x=1.0 y=" + (1.0d / rescaleY) 
                + " z=1.0 width=" + width 
                + " height=" + rescaledHeight
                + " depth=" + depth + " interpolation=Bilinear"
                + " average process create title=impOrthoResized");
        ImagePlus imp2 = IJ.getImage();  // TODO: refactor
        imp2.hide();
        int slices = imp2.getStackSize();
        ImageStack stack = imp2.getStack();
        ImageStack padStack = new ImageStack(width, width, 
                imp.getStackSize());
        for (int s = 1; s <= slices; s++) {
            ImageProcessor ip = stack.getProcessor(s);
            int insertStart = width * (((width - rescaledHeight) / 2) + 1);
            int insertEnd = insertStart + width * rescaledHeight;
            ImageProcessor pip = new FloatProcessor(width, width); // to pad
            float[] pix = (float[])((FloatProcessor)ip).getPixels();
            float[] padpix = new float[width * width];
            for (int i = insertStart; i < insertEnd; i++) {
                padpix[i] = pix[i - insertStart];
            }
            pip.setPixels(padpix);
            padStack.setProcessor(pip, s);
        }
        imp2.setStack(padStack);
        I1l.copyCal(imp, imp2);  // imp2 cal needs to be updated later
        return imp2;
    }
    
    /** Make radial profile plot for each channel at central Z. */
    private ImagePlus makeRadialProfiles(ImagePlus imp) {
        int nc = imp.getNChannels();
        int nz = imp.getNSlices();
        ImagePlus[] profiles = new ImagePlus[nc];
        RadialProfile radialProfiler = new RadialProfile();
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
    
    /** 
     * Overlay resolution rings on each slice of a Fourier ImagePlus. 
     *  NB. raw (non-FFT) imp is required for original calibrations.
     */
    private ImagePlus overlayResRings(ImagePlus Fimp, Calibration cal) {
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
    
    /** main method for manual testing */
    public static void main(String[] args) {
        new ImageJ();
        ImagePlus impTest = TestData.recon;
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
        IJ.selectWindow(impTest.getID());
        IJ.runPlugIn(Rec_FourierPlots.class.getName(), "");
    }
}
