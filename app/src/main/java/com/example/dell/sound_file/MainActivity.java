package com.example.dell.sound_file;

import android.Manifest;
import android.app.AlertDialog;
import android.app.Dialog;
import android.content.Context;
import android.content.DialogInterface;
import android.content.Intent;
import android.content.pm.PackageManager;
import android.location.Location;
import android.location.LocationListener;
import android.location.LocationManager;
import android.media.AudioFormat;
import android.media.AudioRecord;
import android.media.MediaRecorder;
import android.os.Bundle;
import android.os.Environment;
import android.provider.Settings;
import android.support.v4.app.ActivityCompat;
import android.support.v7.app.AppCompatActivity;
import android.util.Log;
import android.view.Menu;
import android.view.MenuInflater;
import android.view.MenuItem;
import android.view.View;
import android.widget.Button;
import android.widget.EditText;
import android.widget.LinearLayout;
import android.widget.TextView;
import android.widget.Toast;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.text.SimpleDateFormat;
import java.util.Date;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

public class MainActivity extends AppCompatActivity {

    private final static int RECORDER_AUDIO_ENCODING = AudioFormat.ENCODING_PCM_16BIT;
    private final static int RECORDER_CHANNELS = AudioFormat.CHANNEL_IN_MONO;
    private final static int RECORDER_SAMPLERATE = 44100;
    private final static int BYTES_PER_ELEMENT = 2;
    private final static int BLOCK_SIZE = AudioRecord.getMinBufferSize(
            RECORDER_SAMPLERATE, RECORDER_CHANNELS, RECORDER_AUDIO_ENCODING)
            / BYTES_PER_ELEMENT;
    private final static int BLOCK_SIZE_FFT = 1764;
    private final static int NUMBER_OF_FFT_PER_SECOND = RECORDER_SAMPLERATE
            / BLOCK_SIZE_FFT;
    private final static double FREQRESOLUTION = ((double) RECORDER_SAMPLERATE)
            / BLOCK_SIZE_FFT;
    LocationListener locationListener;
    LocationManager locationManager;
    MyLocationListener L = new MyLocationListener();
    View mChart;
    LinearLayout layout;
    String fileName;
    File Dir, Dir1;
    FileOutputStream fileOutputStream = null;
    RandomAccessFile randomAccessFile;
    int S_NO = 0;
    Button start_button, stop_button, exit_button, btnYes, btnNo;
    TextView textView;
    String gainString;
    String timeLogString;
    String timeDisplayString;
    String[] THIRD_OCTAVE_LABEL = {"16", "20", "25", "31.5", "40", "50", "63", "80", "100", "125", "160", "200", "250", "315", "400", "500",
            "630", "800", "1000", "1250", "1600", "2000", "2500", "3150", "4000", "5000", "6300", "8000", "10000", "12500", "16000", "20000"};
    // check for level
    String levelToShow;
    // Running Leq
    // final private int finalCountTimeDisplay;
    double linearFftAGlobalRunning = 0;
    // SLM min e max
    double dbFftAGlobalMinTemp = 0;
    double dbFftAGlobalMaxTemp = 0;
    int dbFftAGlobalMinFirst = 0;
    int dbFftAGlobalMaxFirst = 0;
    int kkk = 0; // controllo per bontà leq bande: solo se kkk > 10 misurano bene
    String path1;
    String value = " data";
    String latitude, longitude;
    private int start_flag = 0;
    private DoubleFFT_1D fft = null;
    // Terzi d'ottava
    private float[] THIRD_OCTAVE = {16, 20, 25, 31.5f, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500,
            630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000};
    float[] dbBandMax = new float[THIRD_OCTAVE.length];
    float[] dbBandMin = new float[THIRD_OCTAVE.length];
    private long fftCount = 0;
    private double dbFftAGlobalRunning;
    //
    //check for good leq bands: only if kkk> 10 measure well
    private float[] dbBandRunning = new float[THIRD_OCTAVE.length];
    private float[] linearBandRunning = new float[THIRD_OCTAVE.length];
    // grafico SOnogram    SOnogram chart
//    private float[][] dbHistorySonogram = new float[750][THIRD_OCTAVE.length];
    // variabili finali per time display
    private double dbFftAGlobalMax;
    private double dbFftAGlobalMin;
    private double dbATimeDisplay;
    private float dbFftTimeDisplay[] = new float[BLOCK_SIZE_FFT / 2];
    //    // verifica gain
//    private AutomaticGainControl AGC;
//    private boolean  agcEnable0,agcEnable1,agcEnable2,agcEnable3,agcEnable4,agcEnable5,agcEnable6;
    private float dbFftATimeDisplay[] = new float[BLOCK_SIZE_FFT / 2];
    private float[] dbBandTimeDisplay = new float[THIRD_OCTAVE.length];
    private float[] linearBandTimeDisplay = new float[THIRD_OCTAVE.length];
    private Date dateLogStart;
    // grafico SLMHIstory   graphic SLMHIstory
    private float[] dbAHistoryTimeDisplay = new float[750];
    private float[] dbFftAGlobalRunningHistory = new float[750];
    private int timeLog;
    private String timeLogStringMinSec;
    //    private int timeDisplay;
    private double timeDisplay;
    private AudioRecord recorder;
    private Thread recordingThread = null;
    //private DoubleFFT_1D fft = null;
    private boolean isRecording = false;
    private double filter = 0;
    private double[] weightedA = new double[BLOCK_SIZE_FFT];
    private double actualFreq;
    private float gain;

    private void startNoiseRecording() {
        try {
            gain = Float.parseFloat(gainString);
        } catch (Exception e) {
            gain = 0.0f;
        }
        try {
            timeDisplay = Double.parseDouble(timeDisplayString);
        } catch (Exception e) {
            timeDisplay = 0.5;
        }
        try {
            timeLog = Integer.parseInt(timeLogString);
        } catch (Exception e) {
            timeLog = 1;
        }
        final int finalCountTimeDisplay = (int) (timeDisplay * NUMBER_OF_FFT_PER_SECOND);
        final int finalCountTimeLog = timeLog * NUMBER_OF_FFT_PER_SECOND;

        precalculateWeightedA();

        startRecording(gain, finalCountTimeDisplay, finalCountTimeLog);
    }

    public void setValue(String voice_data, String latitude, String longitude) {
        value = voice_data;
        //
        textView.setText(voice_data);
        S_NO++;
        File Root = Environment.getExternalStorageDirectory();
        try {
            randomAccessFile = new RandomAccessFile(Root + "/" + fileName + "/SOUND.txt", "rw");
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        try {
            randomAccessFile.seek(randomAccessFile.length());
//                randomAccessFile.write(textView1.getText().toString().getBytes(),3,6);
//                randomAccessFile.write(textView2.getText().toString().getBytes(),3,6);
//                randomAccessFile.write(textView3.getText().toString().getBytes(),3,6);
            String store = S_NO + "\t" + fileName + "\t" + voice_data + "\t" + latitude + "\t" + longitude + "\n";
            randomAccessFile.write(store.getBytes());
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void precalculateWeightedA() {
        for (int i = 0; i < BLOCK_SIZE_FFT; i++) {
            double actualFreq = FREQRESOLUTION * i;
            double actualFreqSQ = actualFreq * actualFreq;
            double actualFreqFour = actualFreqSQ * actualFreqSQ;
            double actualFreqEight = actualFreqFour * actualFreqFour;
            double t1 = 20.598997 * 20.598997 + actualFreqSQ;
            t1 = t1 * t1;
            double t2 = 107.65265 * 107.65265 + actualFreqSQ;
            double t3 = 737.86223 * 737.86223 + actualFreqSQ;
            double t4 = 12194.217 * 12194.217 + actualFreqSQ;
            t4 = t4 * t4;

            double weightFormula = (3.5041384e16 * actualFreqEight)
                    / (t1 * t2 * t3 * t4);

            weightedA[i] = weightFormula;
        }
    }

    public void startRecording(final float gain, final int finalCountTimeDisplay, final int finalCountTimeLog) {
        recorder = new AudioRecord(MediaRecorder.AudioSource.VOICE_RECOGNITION,
                RECORDER_SAMPLERATE, RECORDER_CHANNELS,
                RECORDER_AUDIO_ENCODING, BLOCK_SIZE * BYTES_PER_ELEMENT);
        if (recorder.getState() == 1)
            Log.d("our log", "\n" +
                    "The recorder is ready");
        else
            Log.d("our log", "\n" +
                    "The recorder is not ready");

        recorder.startRecording();
        isRecording = true;

        // I create a fft from BLOCK_SIZE_FFT points -> BLOCK_SIZE_FFT / 2 useful bands,
        //each from FREQRESOLUTION Hz
        fft = new DoubleFFT_1D(BLOCK_SIZE_FFT);
        recordingThread = new Thread(new Runnable() {
            public void run() {
                //Raw Data Array (tot: BLOCK_SIZE_FFT * 2 bytes)
                short rawData[] = new short[BLOCK_SIZE_FFT];
                //Unweighed Magnet Array (BLOCK_SIZE_FFT / 2 because it is the number of
                // useful bands)
                final float dbFft[] = new float[BLOCK_SIZE_FFT / 2];
                // Hex core arrays
                final float dbFftA[] = new float[BLOCK_SIZE_FFT / 2];
                float normalizedRawData;
                //The fft works with double and complex numbers (re + im in
                // sequence)
                double[] audioDataForFFT = new double[BLOCK_SIZE_FFT * 2];
                //Audience Threshold (20 * 10 ^ (- 6))
                float amplitudeRef = 0.00002f;
                //third octave
                final float[] dbBand = new float[THIRD_OCTAVE.length];
                final float[] linearBand = new float[THIRD_OCTAVE.length];
                final float[] linearBandCount = new float[THIRD_OCTAVE.length];
                int n = 3;
//                float summingLinearBand = 0f;
//                int controllo_frequenze = 0;
//                int controllo_frequenze_1 = 0;

                // Variables for Medium Time Calculator
                int indexTimeDisplay = 1;
                double linearATimeDisplay = 0;


                // Variables for Time Calculation Time Log
                int indexTimeLog = 0;
                double linearTimeLog = 0;
                double linearATimeLog = 0;
                final float[] linearBandTimeLog = new float[THIRD_OCTAVE.length];

                final float linearFftTimeDisplay[] = new float[BLOCK_SIZE_FFT / 2];
                final float linearFftATimeDisplay[] = new float[BLOCK_SIZE_FFT / 2];

                int initial_delay = 0;
                while (isRecording) {

                    //I read the data
                    recorder.read(rawData, 0, BLOCK_SIZE_FFT);
                    // an initial delay was introduced because at activation there were very high levels of running leq (> 100 dB) and low low (10 dB) due perhaps to the initial startup of the peripheral
                    initial_delay++;
                    if (initial_delay > 20) {
                        for (int i = 0, j = 0; i < BLOCK_SIZE_FFT; i++, j += 2) {
                            // Range [-1,1]
                            normalizedRawData = (float) rawData[i]
                                    / (float) Short.MAX_VALUE;
                            //filter = ((double) (fastA * normalizedRawData))
                            // + (fastB * filter);
                            filter = normalizedRawData;
                            //Hannings window
                            double x = (2 * Math.PI * i) / (BLOCK_SIZE_FFT - 1);
                            double winValue = (1 - Math.cos(x)) * 0.5d;
                            // Real Part
                            audioDataForFFT[j] = filter * winValue;
                            // Imaginary part
                            audioDataForFFT[j + 1] = 0.0;
                        }
                        // FFT
                        fft.complexForward(audioDataForFFT);
                        // Magsum not weighed
                        double linearFftGlobal = 0;
                        //Magsum weighed
                        double linearFftAGlobal = 0;
                        //index for a third octave
                        int k = 0;
                        for (int ki = 0; ki < THIRD_OCTAVE.length; ki++) {
                            linearBandCount[ki] = 0;
                            linearBand[ki] = 0;
                            dbBand[ki] = 0;
                        }

                        //I read up to BLOCK_SIZE_FFT / 2 because in tot I have BLOCK_SIZE_FFT / 2
                        //useful bands
                        for (int i = 0, j = 0; i < BLOCK_SIZE_FFT / 2; i++, j += 2) {

                            double re = audioDataForFFT[j];
                            double im = audioDataForFFT[j + 1];

                            // Magnitude
                            double mag = Math.sqrt((re * re) + (im * im));

                            // Weighted A
                            //to understand: for i = 0 is an invalid value (perhaps less infinite), but does it make sense?
                            //this is found in the graph:
                            // for i = 0 the unweighted has a value, while weighing has not ...
                            double weightFormula = weightedA[i];
                            dbFft[i] = (float) (10 * Math.log10(mag * mag
                                    / amplitudeRef))
                                    + gain;
                            dbFftA[i] = (float) (10 * Math.log10(mag * mag
                                    * weightFormula
                                    / amplitudeRef))
                                    + gain;
                            linearFftGlobal += Math.pow(10, dbFft[i] / 10f);
                            linearFftAGlobal += Math.pow(10, dbFftA[i] / 10f);
                            float linearFft = (float) Math.pow(10, dbFft[i] / 10f);
                            if ((0 <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 17.8f)) {
                                linearBandCount[0] += 1;
                                linearBand[0] += linearFft;
                                dbBand[0] = (float) (10 * Math.log10(linearBand[0]));
                            }
                            if ((17.8f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 22.4f)) {
                                linearBandCount[1] += 1;
                                linearBand[1] += linearFft;
                                dbBand[1] = (float) (10 * Math.log10(linearBand[1]));
                            }
                            if ((22.4f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 28.2f)) {
                                linearBandCount[2] += 1;
                                linearBand[2] += linearFft;
                                dbBand[2] = (float) (10 * Math.log10(linearBand[2]));
                            }
                            if ((28.2f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 35.5f)) {
                                linearBandCount[3] += 1;
                                linearBand[3] += linearFft;
                                dbBand[3] = (float) (10 * Math.log10(linearBand[3]));
                            }
                            if ((35.5f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 44.7f)) {
                                linearBandCount[4] += 1;
                                linearBand[4] += linearFft;
                                dbBand[4] = (float) (10 * Math.log10(linearBand[4]));
                            }
                            if ((44.7f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 56.2f)) {
                                linearBandCount[5] += 1;
                                linearBand[5] += linearFft;
                                dbBand[5] = (float) (10 * Math.log10(linearBand[5]));
                            }
                            if ((56.2f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 70.8f)) {
                                linearBandCount[6] += 1;
                                linearBand[6] += linearFft;
                                dbBand[6] = (float) (10 * Math.log10(linearBand[6]));
                            }
                            if ((70.8f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 89.1f)) {
                                linearBandCount[7] += 1;
                                linearBand[7] += linearFft;
                                dbBand[7] = (float) (10 * Math.log10(linearBand[7]));
                            }
                            if ((89.1f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 112f)) {
                                linearBandCount[8] += 1;
                                linearBand[8] += linearFft;
                                dbBand[8] = (float) (10 * Math.log10(linearBand[8]));
                            }
                            if ((112f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 141f)) {
                                linearBandCount[9] += 1;
                                linearBand[9] += linearFft;
                                dbBand[9] = (float) (10 * Math.log10(linearBand[9]));
                            }
                            if ((141f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 178f)) {
                                linearBandCount[10] += 1;
                                linearBand[10] += linearFft;
                                dbBand[10] = (float) (10 * Math.log10(linearBand[10]));
                            }
                            if ((178f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 224f)) {
                                linearBandCount[11] += 1;
                                linearBand[11] += linearFft;
                                dbBand[11] = (float) (10 * Math.log10(linearBand[11]));
                            }
                            if ((224f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 282f)) {
                                linearBandCount[12] += 1;
                                linearBand[12] += linearFft;
                                dbBand[12] = (float) (10 * Math.log10(linearBand[12]));
                            }
                            if ((282f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 355f)) {
                                linearBandCount[13] += 1;
                                linearBand[13] += linearFft;
                                dbBand[13] = (float) (10 * Math.log10(linearBand[13]));
                            }
                            if ((355f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 447f)) {
                                linearBandCount[14] += 1;
                                linearBand[14] += linearFft;
                                dbBand[14] = (float) (10 * Math.log10(linearBand[14]));
                            }
                            if ((447f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 562f)) {
                                linearBandCount[15] += 1;
                                linearBand[15] += linearFft;
                                dbBand[15] = (float) (10 * Math.log10(linearBand[15]));
                            }
                            if ((562f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 708f)) {
                                linearBandCount[16] += 1;
                                linearBand[16] += linearFft;
                                dbBand[16] = (float) (10 * Math.log10(linearBand[16]));
                            }
                            if ((708f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 891f)) {
                                linearBandCount[17] += 1;
                                linearBand[17] += linearFft;
                                dbBand[17] = (float) (10 * Math.log10(linearBand[17]));
                            }
                            if ((891f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 1122f)) {
                                linearBandCount[18] += 1;
                                linearBand[18] += linearFft;
                                dbBand[18] = (float) (10 * Math.log10(linearBand[18]));
                            }
                            if ((1122f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 1413f)) {
                                linearBandCount[19] += 1;
                                linearBand[19] += linearFft;
                                dbBand[19] = (float) (10 * Math.log10(linearBand[19]));
                            }
                            if ((1413f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 1778f)) {
                                linearBandCount[20] += 1;
                                linearBand[20] += linearFft;
                                dbBand[20] = (float) (10 * Math.log10(linearBand[20]));
                            }
                            if ((1778f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 2239f)) {
                                linearBandCount[21] += 1;
                                linearBand[21] += linearFft;
                                dbBand[21] = (float) (10 * Math.log10(linearBand[21]));
                            }
                            if ((2239f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 2818f)) {
                                linearBandCount[22] += 1;
                                linearBand[22] += linearFft;
                                dbBand[22] = (float) (10 * Math.log10(linearBand[22]));
                            }
                            if ((2818f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 3548f)) {
                                linearBandCount[23] += 1;
                                linearBand[23] += linearFft;
                                dbBand[23] = (float) (10 * Math.log10(linearBand[23]));
                            }
                            if ((3548f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 4467f)) {
                                linearBandCount[24] += 1;
                                linearBand[24] += linearFft;
                                dbBand[24] = (float) (10 * Math.log10(linearBand[24]));
                            }
                            if ((4467f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 5623f)) {
                                linearBandCount[25] += 1;
                                linearBand[25] += linearFft;
                                dbBand[25] = (float) (10 * Math.log10(linearBand[25]));
                            }
                            if ((5623f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 7079f)) {
                                linearBandCount[26] += 1;
                                linearBand[26] += linearFft;
                                dbBand[26] = (float) (10 * Math.log10(linearBand[26]));
                            }
                            if ((7079f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 8913f)) {
                                linearBandCount[27] += 1;
                                linearBand[27] += linearFft;
                                dbBand[27] = (float) (10 * Math.log10(linearBand[27]));
                            }
                            if ((8913f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 11220f)) {
                                linearBandCount[28] += 1;
                                linearBand[28] += linearFft;
                                dbBand[28] = (float) (10 * Math.log10(linearBand[28]));
                            }
                            if ((11220f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 14130f)) {
                                linearBandCount[29] += 1;
                                linearBand[29] += linearFft;
                                dbBand[29] = (float) (10 * Math.log10(linearBand[29]));
                            }
                            if ((14130f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 17780f)) {
                                linearBandCount[30] += 1;
                                linearBand[30] += linearFft;
                                dbBand[30] = (float) (10 * Math.log10(linearBand[30]));
                            }
                            if ((17780f <= i * FREQRESOLUTION) && (i * FREQRESOLUTION < 22390f)) {
                                linearBandCount[31] += 1;
                                linearBand[31] += linearFft;
                                dbBand[31] = (float) (10 * Math.log10(linearBand[31]));
                            }

                        }
                        final double dbFftAGlobal = 10 * Math.log10(linearFftAGlobal);
                        //
                        //   calculation of min and max global FFT weighted value A
                        if (dbFftAGlobal > 0) {
                            if (dbFftAGlobalMinFirst == 0) {
                                dbFftAGlobalMinTemp = dbFftAGlobal;
                                dbFftAGlobalMinFirst = 1;
                            } else {
                                if (dbFftAGlobalMinTemp > dbFftAGlobal) {
                                    dbFftAGlobalMinTemp = dbFftAGlobal;
                                }
                            }
                            if (dbFftAGlobalMaxFirst == 0) {
                                dbFftAGlobalMaxTemp = dbFftAGlobal;
                                dbFftAGlobalMaxFirst = 1;
                            } else {
                                if (dbFftAGlobalMaxTemp < dbFftAGlobal) {
                                    dbFftAGlobalMaxTemp = dbFftAGlobal;
                                }
                            }
                        }
                        dbFftAGlobalMin = dbFftAGlobalMinTemp;
                        dbFftAGlobalMax = dbFftAGlobalMaxTemp;


                        // Running Leq
                        fftCount++;
                        linearFftAGlobalRunning += linearFftAGlobal;
                        dbFftAGlobalRunning = 10 * Math.log10(linearFftAGlobalRunning / fftCount);

                        for (int ki = 0; ki < THIRD_OCTAVE.length; ki++) {
                            linearBandRunning[ki] += linearBand[ki];
                            dbBandRunning[ki] = 10 * (float) Math.log10(linearBandRunning[ki] / fftCount);
                        }

                        // min and max calculation for unweighted dbBand
                        //  I define minimum and maximum bandwidths
                        for (int kk = 0; kk < dbBand.length; kk++) {
                            if (dbBandMax[kk] < dbBand[kk]) {
                                dbBandMax[kk] = dbBand[kk];
                            }
                            if (kkk >= 10) { // check for good leq bands: only if kkk> 10 measure well
                                if (dbBandMin[kk] == 0f) {
                                    if (dbBand[kk] > 0) {
                                        dbBandMin[kk] = dbBand[kk];
                                    }
                                } else if (dbBandMin[kk] > dbBand[kk]) {
                                    dbBandMin[kk] = dbBand[kk];
                                }
                            }
                        }
                        kkk++;


                        // LAeqTimeDisplay
                        //
                        // Calculate Average for Time Display and update charts
                        linearATimeDisplay += linearFftAGlobal;
                        for (int i = 0; i < THIRD_OCTAVE.length; i++) {
                            linearBandTimeDisplay[i] += linearBand[i];
                        }

                        for (int i = 0; i < dbFftTimeDisplay.length; i++) {
                            linearFftTimeDisplay[i] += Math.pow(10, dbFft[i] / 10f);
                            linearFftATimeDisplay[i] += Math.pow(10, dbFftA[i] / 10f);
                        }

                        if (indexTimeDisplay < finalCountTimeDisplay) {
                            indexTimeDisplay++;
                        } else {
                            //
                            //  update data for plot of octave thirds
                            for (int i = 0; i < THIRD_OCTAVE.length; i++) {
                                dbBandTimeDisplay[i] = 10 * (float) Math.log10(linearBandTimeDisplay[i] / finalCountTimeDisplay);
                                linearBandTimeDisplay[i] = 0;
                            }

                            // FFT plot
                            for (int i = 0; i < dbFftTimeDisplay.length; i++) {
                                dbFftTimeDisplay[i] = 10 * (float) Math.log10(linearFftTimeDisplay[i] / finalCountTimeDisplay);
                                dbFftATimeDisplay[i] = 10 * (float) Math.log10(linearFftATimeDisplay[i] / finalCountTimeDisplay);
                                linearFftTimeDisplay[i] = 0;
                                linearFftATimeDisplay[i] = 0;
                            }

                            // timeDisplay data and notification icons
                            dbATimeDisplay = 10 * Math.log10(linearATimeDisplay / finalCountTimeDisplay);
                            indexTimeDisplay = 1;
                            linearATimeDisplay = 0;


                            /*Thread thread = new Thread() {
                                @Override
                                public void run() {
                                    Log.v("NOISE", String.valueOf(dbATimeDisplay));
                                  //  setValue(String.valueOf(dbATimeDisplay));
                                    textView.setText(String.valueOf(dbATimeDisplay));


                                }
                            };
                            thread.start();*/
                            runOnUiThread(new Runnable() {
                                @Override
                                public void run() {
                                    Log.v("NOISE", String.valueOf(dbATimeDisplay));
                                    if (start_flag == 1 && latitude != null && longitude != null) {
                                        setValue(String.valueOf(dbATimeDisplay), latitude, longitude);
                                        //textView.setText(String.valueOf(dbATimeDisplay));
                                    }

                                }
                            });
                        }
                    }
                } // while
            }
        }, "AudioRecorder Thread");
        recordingThread.start();
    }

    private void openChart() {


    }

    public void stopRecording() {
        // stops the recording activity
        if (recorder != null) {
            isRecording = false;
            try {
                recordingThread.join();
                //fos.close();
            } catch (Exception e) {
                Log.d("our log",
                        "The Main Thread can not wait for the secondary audio thread to close");
            }
            recorder.stop();
            recorder.release();
            recorder = null;
            recordingThread = null;
        }
    }

    @Override
    protected void onCreate(Bundle savedInstanceState) {
        super.onCreate(savedInstanceState);
        setContentView(R.layout.activity_main);

        start_button = findViewById(R.id.start);
        stop_button = findViewById(R.id.stop);
        exit_button = findViewById(R.id.exit);
        textView = findViewById(R.id.textView);
        stop_button.setVisibility(View.VISIBLE);

        locationListener = new MyLocationListener();
        locationManager = (LocationManager) getSystemService(Context.LOCATION_SERVICE);
        if (ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_FINE_LOCATION) != PackageManager.PERMISSION_GRANTED && ActivityCompat.checkSelfPermission(this, Manifest.permission.ACCESS_COARSE_LOCATION) != PackageManager.PERMISSION_GRANTED) {
            // TODO: Consider calling
            //    ActivityCompat#requestPermissions
            // here to request the missing permissions, and then overriding
            //   public void onRequestPermissionsResult(int requestCode, String[] permissions,
            //                                          int[] grantResults)
            // to handle the case where the user grants the permission. See the documentation
            // for ActivityCompat#requestPermissions for more details.
            return;
        }
        locationManager.requestLocationUpdates(LocationManager.NETWORK_PROVIDER, 250, 1, locationListener);
        //  Location location = locationManager.getLastKnownLocation(LocationManager.GPS_PROVIDER);
        //File Root = Environment.getExternalStorageDirectory();
        layout = findViewById(R.id.layout);
//        try {
//            randomAccessFile = new RandomAccessFile(Root+"/myfile.txt","rw");
//        } catch (FileNotFoundException e) {
//            e.printStackTrace();
//        }

        SimpleDateFormat formatter = new SimpleDateFormat("yyyy_MM_dd_HH:mm:ss");
        Date now = new Date();
        fileName = formatter.format(now);
        File file = new File(Environment.getExternalStorageDirectory() + "/" + fileName);
        file.mkdir();

        //Dir = new File(Root.getAbsoluteFile() + "/MyFile");
        //File Root2 = Environment.getExternalStorageDirectory();
        //Dir1 = new File(Root2.getAbsoluteFile() + "/MyFile2");
        textView.post(new Runnable() {
            public void run() {
                textView.setText(value);
            }
        });
        start_button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                start_button.setVisibility(View.VISIBLE);
                stop_button.setVisibility(View.VISIBLE);
                startNoiseRecording();
                start_flag = 1;

            }
        });
        stop_button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View v) {
                stop_button.setVisibility(View.VISIBLE);
                start_button.setVisibility(View.VISIBLE);
                stopRecording();
                start_flag = 0;
            }
        });

        exit_button.setOnClickListener(new View.OnClickListener() {
            @Override
            public void onClick(View view) {
                AlertDialog.Builder builder = new AlertDialog.Builder(MainActivity.this);
                builder.setTitle("Dialog");
                builder.setIcon(R.mipmap.ic_launcher);
                builder.setMessage("Do You Want To Close").setCancelable(true).setPositiveButton("OK", new DialogInterface.OnClickListener() {
                    @Override
                    public void onClick(DialogInterface dialogInterface, int i) {
                        finish();
                        Toast.makeText(MainActivity.this, "You Have Pressed Ok", Toast.LENGTH_SHORT).show();
                    }
                }).setNegativeButton("Cancle", new DialogInterface.OnClickListener() {
                    @Override
                    public void onClick(DialogInterface dialogInterface, int i) {
                        dialogInterface.dismiss();
                    }
                });
                AlertDialog alert = builder.create();
                alert.show();
            }
        });

    }

    @Override
    public boolean onCreateOptionsMenu(Menu menu) {
        MenuInflater inflater = getMenuInflater();
        inflater.inflate(R.menu.filemenu, menu);
        return true;
    }

    @Override
    public boolean onOptionsItemSelected(MenuItem item) {
        switch (item.getItemId()) {
            case R.id.donationcafe:
                Toast.makeText(this, "Donation Cafe", Toast.LENGTH_SHORT).show();
                return true;
            case R.id.Settings:
                Toast.makeText(this, "Setting", Toast.LENGTH_SHORT).show();
                return true;
            case R.id.startlogging:
                Toast.makeText(this, "startLogging", Toast.LENGTH_SHORT).show();
                final Dialog dialog = new Dialog(this);
                dialog.setContentView(R.layout.custom_dialog);
                dialog.setTitle("Log File");
                EditText Edit = dialog.findViewById(R.id.editText1);
                Button canle = dialog.findViewById(R.id.dial_cancle);
                Button ok = dialog.findViewById(R.id.dial_ok);
                // if button is clicked, close the custom dialog
                canle.setOnClickListener(new View.OnClickListener() {
                    @Override
                    public void onClick(View view) {
                        dialog.dismiss();
                    }
                });
            case R.id.Log_List:
                Toast.makeText(this, "Log List", Toast.LENGTH_SHORT).show();
            default:
                return false;
        }
    }

    private class MyLocationListener implements LocationListener {
        public void onLocationChanged(Location location) {
            latitude = "" + location.getLatitude();
            longitude = "" + location.getLongitude();
            if (start_flag == 1) {
                //setValue(msg);
                textView.setText(latitude + " " + longitude);
                //Toast.makeText(getBaseContext(), "location lis has been working nnow" + msg, Toast.LENGTH_SHORT).show();
            }
        }


        @Override
        public void onStatusChanged(String s, int i, Bundle bundle) {

        }

        @Override
        public void onProviderEnabled(String s) {
            Toast.makeText(getBaseContext(), "GPS is turned ON!!", Toast.LENGTH_SHORT).show();
        }

        @Override
        public void onProviderDisabled(String s) {
            Intent intent = new Intent(Settings.ACTION_LOCATION_SOURCE_SETTINGS);
            startActivity(intent);
            Toast.makeText(getBaseContext(), "GPS is turned off!!", Toast.LENGTH_SHORT).show();
        }

    }
}
