/*
Copyright (c) 2007-2010 Michael Specht

This file is part of qTrace.

qTrace is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

qTrace is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with qTrace.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Quantifier.h"
#include <QtCore>
#include <QtSvg>
#include <math.h> 
#include <limits>
#include "Tango.h"

#define ABSENCE_MASS_ACCURACY_FACTOR 2.0


k_Quantifier::k_Quantifier(QStringList& ak_Arguments, QSet<r_Parameter::Enumeration> ak_Parameters, QString as_ProgramName, QString as_AdditionalArguments)
	: k_QuantifierBase(ak_Arguments, ak_Parameters, as_ProgramName, as_AdditionalArguments)
{
    parseArguments(ak_Arguments);
}


k_Quantifier::~k_Quantifier()
{
}


void k_Quantifier::run()
{
    if (mk_Peptides.empty() || mk_SpectraFiles.empty())
        printUsageAndExit();
    
    removeNonPeptides(mk_Peptides);
    
    printf("Calculating target peaks for %d peptides...", mk_Peptides.size());
    fflush(stdout);
    
    // create target peaks
    if (mb_UseIsotopeEnvelopes)
    {
        foreach (QString ls_Peptide, mk_Peptides)
        {
            QHash<QString, int> lk_Composition = compositionForPeptide(ls_Peptide);
            double ld_BaseMass = mk_IsotopeEnvelope.massForComposition(lk_Composition);
            
            tk_IsotopeEnvelope lk_Envelope_[2];
            lk_Envelope_[0] = lightEnvelopeForPeptide(ls_Peptide);
            lk_Envelope_[1] = heavyEnvelopeForPeptide(ls_Peptide);
            
            tk_IsotopeEnvelope lk_EnvelopeNormalized_[2];
            for (int i = 0; i < 2; ++i)
                lk_EnvelopeNormalized_[i] = k_IsotopeEnvelope::normalize(lk_Envelope_[i]);
            
            for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
            {
                for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
                {
//                     printf("%s (E%d, %d+)\n", ls_Peptide.toStdString().c_str(), li_Envelope, li_Charge);
                    QString ls_PeptideChargeWeight = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);
                    mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight] = r_EnvelopePeaks();
                    
                    // add forbidden peak if light envelope
                    if (mb_CheckForbiddenPeak && (li_Envelope == 0))
                    {
                        double ld_Mz = (ld_BaseMass -md_HydrogenMass + md_HydrogenMass * li_Charge) / li_Charge;
                        // oy, it's the forbidden peak!
                        mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].mk_ForbiddenIds.insert(mk_Targets.size());
                        mk_Targets.insert(ld_Mz, mk_Targets.size());
/*                        double error = ld_Mz * (5.0 / 1000000.0);
                        double min = ld_Mz - error;
                        double max = ld_Mz + error;
                        printf("FOR %3d: %9.4f (%9.4f - %9.4f)\n", mk_Targets.size() - 1, ld_Mz, min, max);*/
                    }

                    // add peaks from isotope envelope
                    for (int li_PeakIndex = 0; li_PeakIndex < lk_Envelope_[li_Envelope].size(); ++li_PeakIndex)
                    {
                        double ld_Abundance = lk_Envelope_[li_Envelope][li_PeakIndex].first;
                        double ld_NormalizedAbundance = lk_EnvelopeNormalized_[li_Envelope][li_PeakIndex].first;
                        double ld_MassShift = lk_Envelope_[li_Envelope][li_PeakIndex].second;
                        double ld_Mz = (ld_BaseMass + ld_MassShift + md_HydrogenMass * li_Charge) / li_Charge;
                        if (ld_NormalizedAbundance >= md_RequireAbundance)
                        {
                            // oy, it's a required peak!
                            mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].mk_RequiredIds.insert(mk_Targets.size());
                            mk_TargetIntensity[mk_Targets.size()] = ld_Abundance;
                            mk_Targets.insert(ld_Mz, mk_Targets.size());
/*                            double error = ld_Mz * (5.0 / 1000000.0);
                            double min = ld_Mz - error;
                            double max = ld_Mz + error;
                            printf("REQ %3d: %9.4f (%9.4f - %9.4f)\n", mk_Targets.size() - 1, ld_Mz, min, max);*/
                        }
                        else if (ld_NormalizedAbundance >= md_ConsiderAbundance)
                        {
                            // oy, it's a peak to be considered!
                            mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].mk_ConsideredIds.insert(mk_Targets.size());
                            mk_TargetIntensity[mk_Targets.size()] = ld_Abundance;
                            mk_Targets.insert(ld_Mz, mk_Targets.size());
/*                            double error = ld_Mz * (5.0 / 1000000.0);
                            double min = ld_Mz - error;
                            double max = ld_Mz + error;
                            printf("CON %3d: %9.4f (%9.4f - %9.4f)\n", mk_Targets.size() - 1, ld_Mz, min, max);*/
                        }
                    }
                }
            }
        }
    }
    
    printf("done.\n");
    fflush(stdout);
    
/*    for (int i = 0; i < mk_Targets.size(); ++i)
        printf("%3d %9.4f %3d\n", i, mk_Targets.keys()[i], mk_Targets.values()[i]);*/
    
    if (mk_pCsvStream)
        *(mk_pCsvStream.get_Pointer()) << "Filename,Scan id,Peptide,Charge,Amount light,Amount heavy,Retention time\n";
    
    // parse all spectra files
    foreach (QString ls_Path, mk_SpectraFiles)
    {
        ms_CurrentSpectraFile = QFileInfo(ls_Path).baseName();

        this->parseFile(ls_Path);
        
        if (!mb_Quiet)
            printf(" done.\n");
    }
}


void k_Quantifier::handleScan(r_Scan& ar_Scan, bool& ab_Continue)
{
    if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
    {
        printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
        return;
    }

	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum, md_MinSnr);
    QList<double> lk_AllPeakMz;
    foreach (r_Peak lr_Peak, lk_AllPeaks)
//     {
        lk_AllPeakMz.append(lr_Peak.md_PeakMz);
//         printf("%9.4f ", lr_Peak.md_PeakMz);
//     }
//     printf("\n");
    
    QHash<int, int> lk_PreMatches = matchTargetsToPeaks(lk_AllPeakMz, mk_Targets, md_MassAccuracy);
    QHash<int, int> lk_Matches;
    foreach (int li_Key, lk_PreMatches.keys())
    {
/*        printf("%d => %d (%d => %d)\n", 
               li_Key, lk_PreMatches[li_Key], 
               mk_Targets.values()[li_Key], lk_PreMatches[li_Key]);*/
        lk_Matches[mk_Targets.values()[li_Key]] = lk_PreMatches[li_Key];
    }
    
    QSet<int> lk_MatchedTargetIds = lk_Matches.keys().toSet();
    
    foreach (QString ls_Peptide, mk_Peptides)
    {
        for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
        {
            if (mb_UseIsotopeEnvelopes)
            {
                QList<tk_DoublePair> lk_Pairs_[2] = {QList<tk_DoublePair>(), QList<tk_DoublePair>()};
                double ld_FitFactor_[2] = {0.0, 0.0};
                double ld_FitError_[2] = {0.0, 0.0};
                bool lb_Good_[2] = {false, false};
                
                for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
                {
                    QString ls_PeptideChargeWeightKey = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);

                    // don't quantify this at all if the forbidden peak was there
                    if ((li_Envelope == 0) && (mb_CheckForbiddenPeak && (!(mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey].mk_ForbiddenIds & lk_MatchedTargetIds).empty())))
                    {
                        lb_Good_[0] = lb_Good_[1] = false;
                        break;
                    }
                    
                    // skip this if not all required light and heavy peaks are there
                    if ((lk_MatchedTargetIds & mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey].mk_RequiredIds).size() == mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey].mk_RequiredIds.size())
                    {
                        QString ls_PeptideChargeWeightKey = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);

                        QSet<int> lk_AllTargets = 
                            mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey].mk_RequiredIds | 
                            mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey].mk_ConsideredIds;
                        foreach (int li_Id, lk_AllTargets)
                        {
                            if (lk_MatchedTargetIds.contains(li_Id))
                            {
                                double ld_TargetIntensity = mk_TargetIntensity[li_Id];
                                double ld_PeakIntensity = lk_AllPeaks[lk_Matches[li_Id]].md_PeakIntensity;
                                lk_Pairs_[li_Envelope] << tk_DoublePair(ld_TargetIntensity, ld_PeakIntensity);
                            }
                        }
                        leastSquaresFit(lk_Pairs_[li_Envelope], &ld_FitFactor_[li_Envelope], &ld_FitError_[li_Envelope]);
                        if (ld_FitError_[li_Envelope] > md_MaxFitError)
                        {
                            lb_Good_[li_Envelope] = false;
                        }
                        else
                            lb_Good_[li_Envelope] = true;
                    }
                }
                if (lb_Good_[0] || lb_Good_[1])
                {
                    if (mk_pCsvStream)
                    {
                        *(mk_pCsvStream.get_Pointer()) << QString("\"%1\",\"%2\",%3,%4,%5,%6,%7\n")
                            .arg(ms_CurrentSpectraFile)
                            .arg(ar_Scan.ms_Id)
                            .arg(ls_Peptide)
                            .arg(li_Charge)
                            .arg(ld_FitFactor_[0], 1, 'f', 4)
                            .arg(ld_FitFactor_[1], 1, 'f', 4)
                            .arg(ar_Scan.md_RetentionTime);
                    }
                }
            }
        }
    }
}


void k_Quantifier::parseArguments(QStringList& ak_Arguments)
{   
    int li_Mode = 0;
    while (!ak_Arguments.empty())
    {
        QString ls_Key = ak_Arguments.takeFirst();
        if (ls_Key == "--spectraFiles")
        {
            li_Mode = 1;
            continue;
        } 
        else if (ls_Key == "--peptideFiles")
        {
            li_Mode = 2;
            continue;
        }
        else if (ls_Key == "--peptides")
        {
            li_Mode = 3;
            continue;
        }
        else
        {
            if (li_Mode == 1)
                mk_SpectraFiles.append(ls_Key);
            else if (li_Mode == 2)
            {
                QFile lk_File(ls_Key);
                if (!lk_File.open(QIODevice::ReadOnly))
                {
                    printf("Error: Unable to open %s.\n", ls_Key.toStdString().c_str());
                    exit(1);
                }
                QTextStream lk_Stream(&lk_File);
                while (!lk_Stream.atEnd())
                {
                    QString ls_Peptide = lk_Stream.readLine().trimmed();
                    if (!ls_Peptide.isEmpty())
                        mk_Peptides.insert(ls_Peptide);
                }
                lk_File.close();
            }
            else if (li_Mode == 3)
            {
                QString ls_Peptide = ls_Key.trimmed();
                if (!ls_Peptide.isEmpty())
                    mk_Peptides.insert(ls_Peptide);
            }
        }
    }
}
