/*
Copyright (c) 2007-2008 Michael Specht

This file is part of SimQuant.

SimQuant is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SimQuant is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SimQuant.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <QtCore>
#include <stdio.h>
#include <ptb/FastaReader.h>
#include "Quantifier.h"
#include <ptb/RefPtr.h>
#include "version.h"

#define HYDROGEN_MASS 1.0078250321

struct r_EnvelopePeaks
{
    QSet<int> mk_RequiredIds;
    QSet<int> mk_ConsideredIds;
    QSet<int> mk_ForbiddenIds;
};

k_Quantifier* lk_Quantifier_;
QString ls_Label = "15N";
r_ScanType::Enumeration le_ScanType = r_ScanType::All;
r_AmountEstimation::Enumeration le_AmountEstimation = r_AmountEstimation::Intensity;
int li_MinCharge = 2;
int li_MaxCharge = 3;
double ld_MinSnr = 2.0;
double ld_MassAccuracy = 5.0;
double ld_RequireAbundance = 0.5;
double ld_ConsiderAbundance = 0.05;
double ld_MaxFitError = 0.01;
bool lb_CheckForbiddenPeak = true;
bool lb_PrintStatusMessages = true;
bool lb_LogScale = true;
bool lb_QuickSearch = true;
int li_EfficiencySearchStart = 9000;
int li_EfficiencySearchEnd = 10000;

QHash<QString, int> lk_Composition;
QHash<int, tk_IsotopeEnvelope> lk_IsotopeEnvelopes;
double ld_BaseMass;
tk_IsotopeEnvelope lk_LightIsotopeEnvelope;


void printUsageAndExit()
{
    printf("Usage: qtrace-estimate [options] [spectra file] [peptide] [retention time]\n");
    printf("The spectra file may be mzData, mzXML or mzML, optionally compressed (.gz|.bz2|.zip).\n");
    printf("Options:\n");
    printf("  --searchType [quick|exhaustive] (default: quick)\n");
    printf("  --efficiencyMin [float] (default: 0.9)\n");
    printf("  --efficiencyMax [float] (default: 1.0)\n");
    printf("  --scanType [full|sim|all] (default: all)\n");
    printf("  --minCharge [int] (default: 2)\n");
    printf("  --maxCharge [int] (default: 3)\n");
    printf("  --minSnr [float] (default: 2.0)\n");
    printf("  --massAccuracy (ppm) [float] (default: 5.0)\n");
    printf("      Specify the mass accuracy that should be used for peak matching.\n");
    printf("  --requireAbundance [float] (default: 0.5)\n");
    printf("      Specify which peaks must be present in an isotope envelope.\n");
    printf("  --considerAbundance [float] (default: 0.05)\n");
    printf("      Specify which peaks may additionally be taken into account.\n");
    printf("  --checkForbiddenPeak [flag] (default: yes)\n");
    printf("      Specify whether the forbidden peak is required to be absent.\n");
    printf("  --csvOutput [flag] (default: yes)\n");
    printf("      Enable or disable CSV output.\n");
    printf("  --csvOutputTarget [path] (default: stdout)\n");
    printf("      Redirect CSV output to a file. Enables CSV output.\n");
    printf("  --xhtmlOutputTarget [path] (default: none)\n");
    printf("      Output quantitation events as XHTML-embedded SVG to [path].\n");
    printf("  --logScale [flag] (default: yes)\n");
    printf("      Use logarithmic scale in XHTML spectra.\n");
    printf("  --quiet\n");
    printf("      Don't print status messages.\n");
    printf("  --version\n");
    printf("      Print version and exit.\n");
    exit(1);
}


bool stringToBool(QString& as_String)
{
    as_String = as_String.toLower().trimmed();
    if (as_String == "yes" || as_String == "true" || as_String == "on" || as_String == "enable" || as_String == "enabled")
        return true;
    else if (as_String == "no" || as_String == "false" || as_String == "off" || as_String == "disable" || as_String == "disabled")
        return false;
    else
    {
        printf("Error: unknown boolean value '%s'.\n", as_String.toStdString().c_str());
        exit(1);
    }
};


// returns fit error or a negative value on 'no can do!'
double determineFitError(QList<r_Peak> ak_Peaks, int ai_Efficiency, int ai_Charge)
{
    double ld_Error = -1.0;
    
    // create target isotope envelope if it doesn't exist
    if (!lk_IsotopeEnvelopes.contains(ai_Efficiency))
    {
        tk_ModifiedAbundances lk_Environment;
        double ld_Efficiency = (double)ai_Efficiency / 10000.0;
        lk_Environment["N"] = QList<double>() << (1.0 - ld_Efficiency) << ld_Efficiency;
        k_IsotopeEnvelope lk_IsotopeEnvelope(lk_Environment);
        tk_IsotopeEnvelope lk_Envelope = lk_IsotopeEnvelope.isotopeEnvelopeForComposition(lk_Composition);
        lk_IsotopeEnvelopes[ai_Efficiency] = lk_Envelope;
    }
    
    r_EnvelopePeaks lr_LightPeptide;
    r_EnvelopePeaks lr_HeavyPeptide;
    
    tk_IsotopeEnvelope& lk_Envelope = lk_LightIsotopeEnvelope;
    tk_IsotopeEnvelope lk_NormalizedEnvelope = k_IsotopeEnvelope::normalize(lk_Envelope);
    
    QMultiMap<double, int> lk_TargetMasses;
    QHash<int, double> lk_IntensityForId;
    QHash<int, double> lk_NormalizedIntensityForId;
    
    int li_Id = 0;

    // add most prominent isotope envelope peaks (light peptide)
    for (int i = 0; i < lk_Envelope.size(); ++i)
    {
        double ld_NormalizedAbundance = lk_NormalizedEnvelope[i].first;
        if (ld_NormalizedAbundance >= ld_ConsiderAbundance)
        {
            double ld_MassShift = lk_NormalizedEnvelope[i].second;
            double ld_Mz = (ld_BaseMass + ld_MassShift + HYDROGEN_MASS * ai_Charge) / ai_Charge;
            ++li_Id;
            lk_TargetMasses.insert(ld_Mz, li_Id);
            lk_IntensityForId[li_Id] = lk_Envelope[i].first;
            lk_NormalizedIntensityForId[li_Id] = lk_NormalizedEnvelope[i].first;
            if (ld_NormalizedAbundance >= ld_RequireAbundance)
                lr_LightPeptide.mk_RequiredIds.insert(li_Id);
            else
                lr_LightPeptide.mk_ConsideredIds.insert(li_Id);
        }
    }
    
    // add forbidden peak (light peptide)
    double ld_Mz = (ld_BaseMass - HYDROGEN_MASS + HYDROGEN_MASS * ai_Charge) / ai_Charge;
    ++li_Id;
    lk_TargetMasses.insert(ld_Mz, li_Id);
    lr_LightPeptide.mk_ForbiddenIds.insert(li_Id);
    
    lk_Envelope = lk_IsotopeEnvelopes[ai_Efficiency];
    lk_NormalizedEnvelope = k_IsotopeEnvelope::normalize(lk_Envelope);

    // add most prominent isotope envelope peaks (heavy peptide)
    for (int i = 0; i < lk_Envelope.size(); ++i)
    {
        double ld_NormalizedAbundance = lk_NormalizedEnvelope[i].first;
        if (ld_NormalizedAbundance >= ld_ConsiderAbundance)
        {
            double ld_MassShift = lk_NormalizedEnvelope[i].second;
            double ld_Mz = (ld_BaseMass + ld_MassShift + HYDROGEN_MASS * ai_Charge) / ai_Charge;
            ++li_Id;
            lk_TargetMasses.insert(ld_Mz, li_Id);
            lk_IntensityForId[li_Id] = lk_Envelope[i].first;
            lk_NormalizedIntensityForId[li_Id] = lk_NormalizedEnvelope[i].first;
            if (ld_NormalizedAbundance >= ld_RequireAbundance)
                lr_HeavyPeptide.mk_RequiredIds.insert(li_Id);
            else
                lr_HeavyPeptide.mk_ConsideredIds.insert(li_Id);
        }
    }
    
    // build peak m/z list
    QList<double> lk_PeakMz;
    foreach (r_Peak lr_Peak, ak_Peaks)
        lk_PeakMz.append(lr_Peak.md_PeakMz);
    QList<double> lk_TargetMzSorted = lk_TargetMasses.keys();
    QList<int> lk_TargetMzIdSorted = lk_TargetMasses.values();
    QHash<int, int> lk_PeakForTarget = lk_Quantifier_->matchTargetsToPeaks(lk_PeakMz, lk_TargetMzSorted, ld_MassAccuracy);
    
    QHash<int, int> lk_MatchedLightRequiredPeaks = lk_Quantifier_->extractMatches(lr_LightPeptide.mk_RequiredIds, lk_TargetMzIdSorted, lk_PeakForTarget);
    QHash<int, int> lk_MatchedLightConsideredPeaks = lk_Quantifier_->extractMatches(lr_LightPeptide.mk_ConsideredIds, lk_TargetMzIdSorted, lk_PeakForTarget);
    QHash<int, int> lk_MatchedLightForbiddenPeaks = lk_Quantifier_->extractMatches(lr_LightPeptide.mk_ForbiddenIds, lk_TargetMzIdSorted, lk_PeakForTarget);
    QHash<int, int> lk_MatchedHeavyRequiredPeaks = lk_Quantifier_->extractMatches(lr_HeavyPeptide.mk_RequiredIds, lk_TargetMzIdSorted, lk_PeakForTarget);
    QHash<int, int> lk_MatchedHeavyConsideredPeaks = lk_Quantifier_->extractMatches(lr_HeavyPeptide.mk_ConsideredIds, lk_TargetMzIdSorted, lk_PeakForTarget);
    
    /*
    printf("Matched %d of %d peaks (L %d/%d %d/%d F %d/1 H %d/%d %d/%d).\n", 
           lk_PeakForTarget.size(), lk_TargetMzSorted.size(),
           lk_MatchedLightRequiredPeaks.size(),
           lr_LightPeptide.mk_RequiredIds.size(),
           lk_MatchedLightConsideredPeaks.size(),
           lr_LightPeptide.mk_ConsideredIds.size(),
           lk_MatchedLightForbiddenPeaks.size(),
           lk_MatchedHeavyRequiredPeaks.size(),
           lr_HeavyPeptide.mk_RequiredIds.size(),
           lk_MatchedHeavyConsideredPeaks.size(),
           lr_HeavyPeptide.mk_ConsideredIds.size());
           */

    if (lb_CheckForbiddenPeak)
    {
        if (!lk_MatchedLightForbiddenPeaks.empty())
            return -1.0;
    }
    
    if (lk_MatchedLightRequiredPeaks.size() == lr_LightPeptide.mk_RequiredIds.size())
    {
        if (lk_MatchedHeavyRequiredPeaks.size() == lr_HeavyPeptide.mk_RequiredIds.size())
        {
            // all required light and heavy peaks are there, this must be the peptide!
            // now do a least squares fit of the observed peaks to the calculated isotope envelope
            QList<tk_DoublePair> lk_Pairs;
            // lk_PairFactor contains a factor from 0 to 1 which denotes the weight of the pair
            // (smaller target peaks are not so important, after all!)
            QList<double> lk_PairFactor;
            double ld_PairFactorSum = 0.0;
            foreach (int li_Id, lr_HeavyPeptide.mk_RequiredIds)
            {
                 double f = pow(lk_NormalizedIntensityForId[li_Id], 0.5);
//                 double f = 1.0;
                lk_PairFactor << f;
                ld_PairFactorSum += f;
                lk_Pairs << tk_DoublePair(lk_IntensityForId[li_Id], ak_Peaks[lk_MatchedHeavyRequiredPeaks[li_Id]].md_PeakIntensity);
            }
            foreach (int li_Id, lr_HeavyPeptide.mk_ConsideredIds)
            {
                if (lk_MatchedHeavyConsideredPeaks.contains(li_Id))
                {
                    double f = pow(lk_NormalizedIntensityForId[li_Id], 0.5);
//                     double f = 1.0;
                    lk_PairFactor << f;
                    ld_PairFactorSum += f;
                    lk_Pairs << tk_DoublePair(lk_IntensityForId[li_Id], ak_Peaks[lk_MatchedHeavyConsideredPeaks[li_Id]].md_PeakIntensity);
                }
            }
            
            double ld_MaxTargetIntensity = 0.0;
            foreach (tk_DoublePair lk_Pair, lk_Pairs)
                ld_MaxTargetIntensity = std::max<double>(ld_MaxTargetIntensity, lk_Pair.first);
            
//             printf("max target intensity: %1.2f\n", ld_MaxTargetIntensity);

            double ld_Factor = lk_Quantifier_->leastSquaresFit(lk_Pairs);
//             printf("factor: %1.2f\n", ld_Factor);
            ld_Error = 0.0;
            // the error is the sum of squared differences (where the target set is normalized to max. intensity 1)

            for (int i = 0; i < lk_Pairs.size(); ++i)
            {
                tk_DoublePair lk_Pair = lk_Pairs[i];
                double ld_PairFactor = lk_PairFactor[i];
                double ld_TargetHeight = lk_Pair.first;
                double ld_PeakHeight = lk_Pair.second;
                ld_PeakHeight /= ld_Factor;
                // normalize so that max target height is 1.0
                ld_TargetHeight /= ld_MaxTargetIntensity;
                ld_PeakHeight /= ld_MaxTargetIntensity;
//                 printf("%1.2f / %1.2f\n", ld_TargetHeight, ld_PeakHeight);
                ld_Error += pow(ld_PeakHeight - ld_TargetHeight, 2.0) * ld_PairFactor;
            }
            ld_Error /= ld_PairFactorSum;
            ld_Error /= lk_Pairs.size();
        }
    }
    
    return ld_Error;
}


int main(int ai_ArgumentCount, char** ac_Arguments__)
{
    QStringList lk_Arguments;
    for (int i = 1; i < ai_ArgumentCount; ++i)
        lk_Arguments << ac_Arguments__[i];
        
    if (!lk_Arguments.empty() && (lk_Arguments.first() == "--version"))
    {
        printf("qtrace-estimate %s\n", gs_Version.toStdString().c_str());
        exit(0);
    }

    QFile lk_StdOut;
    lk_StdOut.open(stdout, QIODevice::WriteOnly);
    
    RefPtr<QFile> lk_pCsvOutFile;
    RefPtr<QFile> lk_pXhtmlOutFile;
    
    QIODevice* lk_CsvDevice_ = &lk_StdOut;
//     QIODevice* lk_XhtmlDevice_ = &lk_StdOut;
    QIODevice* lk_XhtmlDevice_ = NULL;
    
    // consume options
    int li_Index;
    
    li_Index = lk_Arguments.indexOf("--scanType");
    if (li_Index > -1)
    {
        QString ls_ScanType = lk_Arguments[li_Index + 1];
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        if (ls_ScanType == "sim")
            le_ScanType = r_ScanType::SIM;
        else if (ls_ScanType == "full")
        // why would we say MS1 and MSn here? Because I saw MS1 full scans that
        // were marked as MSn with ms level 1, which probably should be MS1 in the
        // first place, but what can a girl do? It just means don't use SIM.
            le_ScanType = (r_ScanType::Enumeration)(r_ScanType::MS1 | r_ScanType::MSn);
        else if (ls_ScanType == "all")
            le_ScanType = r_ScanType::All;
        else
        {
            printf("Error: unknown scan type %s.\n", ls_ScanType.toStdString().c_str());
            exit(1);
        }
    }
    
    li_Index = lk_Arguments.indexOf("--searchType");
    if (li_Index > -1)
    {
        QString ls_Value = lk_Arguments[li_Index + 1];
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        if (ls_Value == "quick")
            lb_QuickSearch = true;
        else if (ls_Value == "exhaustive")
            lb_QuickSearch = false;
        else
        {
            printf("Error: unknown search type %s.\n", ls_Value.toStdString().c_str());
            exit(1);
        }
    }
    
    li_Index = lk_Arguments.indexOf("--efficiencyMin");
    if (li_Index > -1)
    {
        QString ls_Value = lk_Arguments[li_Index + 1];
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        bool lb_Ok = true;
        double ld_Value = ls_Value.toDouble(&lb_Ok);
        if (!lb_Ok)
        {
            printf("Error: Invalid float: %s.\n", ls_Value.toStdString().c_str());
            exit(1);
        }
        if (ld_Value < 0.0 || ld_Value > 1.0)
        {
            printf("Error: minimum labelling efficiency must be within a range from 0 to 1.\n");
            exit(1);
        }
        li_EfficiencySearchStart = (int)(ld_Value * 10000.0);
    }
    
    li_Index = lk_Arguments.indexOf("--efficiencyMax");
    if (li_Index > -1)
    {
        QString ls_Value = lk_Arguments[li_Index + 1];
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
        bool lb_Ok = true;
        double ld_Value = ls_Value.toDouble(&lb_Ok);
        if (!lb_Ok)
        {
            printf("Error: Invalid float: %s.\n", ls_Value.toStdString().c_str());
            exit(1);
        }
        if (ld_Value < 0.0 || ld_Value > 1.0)
        {
            printf("Error: maximum labelling efficiency must be within a range from 0 to 1.\n");
            exit(1);
        }
        li_EfficiencySearchEnd = (int)(ld_Value * 10000.0);
    }
    
    li_Index = lk_Arguments.indexOf("--minCharge");
    if (li_Index > -1)
    {
        li_MinCharge = QVariant(lk_Arguments[li_Index + 1]).toInt();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--maxCharge");
    if (li_Index > -1)
    {
        li_MaxCharge = QVariant(lk_Arguments[li_Index + 1]).toInt();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--minSnr");
    if (li_Index > -1)
    {
        ld_MinSnr = QVariant(lk_Arguments[li_Index + 1]).toDouble();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--massAccuracy");
    if (li_Index > -1)
    {
        ld_MassAccuracy = QVariant(lk_Arguments[li_Index + 1]).toDouble();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--requireAbundance");
    if (li_Index > -1)
    {
        ld_RequireAbundance = QVariant(lk_Arguments[li_Index + 1]).toDouble();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--considerAbundance");
    if (li_Index > -1)
    {
        ld_ConsiderAbundance = QVariant(lk_Arguments[li_Index + 1]).toDouble();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--maxFitError");
    if (li_Index > -1)
    {
        ld_MaxFitError = QVariant(lk_Arguments[li_Index + 1]).toDouble();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--csvOutput");
    if (li_Index > -1)
    {
        QString ls_Value = lk_Arguments[li_Index + 1];
        if (!stringToBool(ls_Value))
            lk_CsvDevice_ = NULL;
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--csvOutputTarget");
    if (li_Index > -1)
    {
        lk_pCsvOutFile = RefPtr<QFile>(new QFile(lk_Arguments[li_Index + 1]));
        lk_pCsvOutFile->open(QIODevice::WriteOnly);
        lk_CsvDevice_ = lk_pCsvOutFile.get_Pointer();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--xhtmlOutputTarget");
    if (li_Index > -1)
    {
        lk_pXhtmlOutFile = RefPtr<QFile>(new QFile(lk_Arguments[li_Index + 1]));
        lk_pXhtmlOutFile->open(QIODevice::WriteOnly);
        lk_XhtmlDevice_ = lk_pXhtmlOutFile.get_Pointer();
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--checkForbiddenPeak");
    if (li_Index > -1)
    {
        QString ls_Value = lk_Arguments[li_Index + 1];
        lb_CheckForbiddenPeak = stringToBool(ls_Value);
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--logScale");
    if (li_Index > -1)
    {
        QString ls_Value = lk_Arguments[li_Index + 1];
        lb_LogScale = stringToBool(ls_Value);
        lk_Arguments.removeAt(li_Index);
        lk_Arguments.removeAt(li_Index);
    }
    
    li_Index = lk_Arguments.indexOf("--quiet");
    if (li_Index > -1)
    {
        lk_Arguments.removeAt(li_Index);
        lb_PrintStatusMessages = false;
    }
    
    //RefPtr<QIODevice> lk_pTextDevice(new QIODevice(stdout));
    
    QStringList lk_SpectraFiles;
    QHash<QString, double> lk_Peptides;
    
    if (lk_Arguments.empty())
        printUsageAndExit();
    
    lk_SpectraFiles.append(lk_Arguments.takeFirst());
    
    QString ls_Peptide;
    double ld_RetentionTime;
    if (lk_Arguments.size() < 2)
    {
        printf("Error: A peptide and its retention time must be specified.\n");
        exit(1);
    }
    ls_Peptide = lk_Arguments.takeFirst();
    bool lb_Ok = false;
    ld_RetentionTime = lk_Arguments.takeFirst().toDouble(&lb_Ok);
    if (!lb_Ok)
    {
        printf("Error: Invalid retention time specified.\n");
        exit(1);
    }
    
    if (lk_SpectraFiles.empty())
        printUsageAndExit();
        
    lk_Quantifier_ = 
        new k_Quantifier(ls_Label, le_ScanType, le_AmountEstimation,
                     QList<tk_IntPair>() << tk_IntPair(1, 1),
                     li_MinCharge, li_MaxCharge, ld_MinSnr, 
                     ld_MassAccuracy, ld_RequireAbundance, 
                     ld_ConsiderAbundance, ld_MaxFitError, lk_CsvDevice_, 
                     lk_XhtmlDevice_, lb_CheckForbiddenPeak, 
                     lb_PrintStatusMessages, lb_LogScale);
                      
    lk_Composition = lk_Quantifier_->compositionForPeptide(ls_Peptide);
    k_IsotopeEnvelope lk_IsotopeEnvelope;
    lk_LightIsotopeEnvelope = lk_IsotopeEnvelope.isotopeEnvelopeForComposition(lk_Composition);
    ld_BaseMass = lk_IsotopeEnvelope.massForComposition(lk_Composition);
        
    tk_ScanList lk_ScanList = lk_Quantifier_->estimate(lk_SpectraFiles, ls_Peptide, ld_RetentionTime);
    
    QList<double> lk_DeterminedEfficiencies;

    // now check the scan list
    foreach (r_Scan lr_Scan, lk_ScanList)
    {
        QList<r_Peak> lk_Peaks = k_ScanIterator::findAllPeaks(lr_Scan.mr_Spectrum);
        for (int li_Charge = li_MinCharge; li_Charge <= li_MaxCharge; ++li_Charge)
        {
            QMap<int, double> lk_ErrorForEfficiency;
            int li_FinalEfficiency = 0;
            int li_EfficiencyStart = li_EfficiencySearchStart;
            int li_EfficiencyEnd = li_EfficiencySearchEnd;
            int li_EfficiencyStep = 100;
            // if we do the exhaustive search, sample every single point directly
            if (!lb_QuickSearch)
                li_EfficiencyStep = 1;
            while (li_EfficiencyStep > 0)
            {
                QMultiMap<double, int> lk_TestResults;
                for (int li_Efficiency = li_EfficiencyStart; li_Efficiency <= li_EfficiencyEnd; li_Efficiency += li_EfficiencyStep)
                {
                    double ld_Error = -1.0;
                    // re-use cached result if possible
                    if (lk_ErrorForEfficiency.contains(li_Efficiency))
                        ld_Error = lk_ErrorForEfficiency[li_Efficiency];
                    else
                        ld_Error = determineFitError(lk_Peaks, li_Efficiency, li_Charge);
                    if (ld_Error >= 0.0)
                    {
                        lk_TestResults.insert(ld_Error, li_Efficiency);
                        lk_ErrorForEfficiency.insert(li_Efficiency, ld_Error);
                    }
//                      printf("%6.2f%% (%d+) efficiency: %e\n", (double)li_Efficiency / 100.0, li_Charge, ld_Error);
                }
                
                // abort the search if nothing found
                if (lk_TestResults.empty())
                    break;
                
                // dig further into the minimum range
                li_FinalEfficiency = lk_TestResults.values().first();

                li_EfficiencyStart = li_FinalEfficiency;
                li_EfficiencyStart -= li_EfficiencyStep;
                li_EfficiencyStart = std::max<int>(li_EfficiencyStart, 0);
                
                li_EfficiencyEnd = li_FinalEfficiency;
                li_EfficiencyEnd += li_EfficiencyStep;
                li_EfficiencyEnd = std::min<int>(li_EfficiencyEnd, 10000);
                
                li_EfficiencyStep /= 10;
            }
            // the determined labeling efficiency is in li_FinalEfficiency, invalid if 0
            if (li_FinalEfficiency > 0)
            {
                lk_DeterminedEfficiencies.append((double)li_FinalEfficiency / 10000.0);
//                 printf("scan #%s (%d+): %1.2f%%\n", lr_Scan.ms_Id.toStdString().c_str(), li_Charge, (double)li_FinalEfficiency / 100.0);
            }
/*                foreach (int li_Efficiency, lk_ErrorForEfficiency.keys())
                    printf("%1.2f,%f\n", (double)li_Efficiency / 100.0, lk_ErrorForEfficiency[li_Efficiency]);*/
        }
    }
    
    double ld_Mean, ld_SD;
    qSort(lk_DeterminedEfficiencies.begin(), lk_DeterminedEfficiencies.end());
    lk_Quantifier_->calculateMeanAndStandardDeviation(lk_DeterminedEfficiencies, &ld_Mean, &ld_SD);
    printf("%1.2f%% (SD %1.2f%%, RSD %1.4f): %s\n", 
           ld_Mean * 100.0, ld_SD * 100.0, ld_SD / ld_Mean,
           ls_Peptide.toStdString().c_str());
    
    delete lk_Quantifier_;
}
