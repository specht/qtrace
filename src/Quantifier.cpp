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
#include <math.h> 
#include <limits>
#include "Tango.h"


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
    {
        if (mk_Peptides.empty())
            printf("Error: No peptides have been specified.\n");
        if (mk_SpectraFiles.empty())
            printf("Error: No spectra files have been specified.\n");
        exit(1);
    }
    
    removeNonPeptides(mk_Peptides);
    
    printf("Calculating target peaks for %d peptide%s... ", mk_Peptides.size(), mk_Peptides.size() == 1 ? "" : "s");
//     /*DEBUG*/printf("\n");
    fflush(stdout);
    
    QStringList lk_StarAminoAcidOrder = mk_StarAminoAcids.toList();
    
    // create target peaks
    foreach (QString ls_Peptide, mk_Peptides)
    {
        QHash<QString, int> lk_Composition = compositionForPeptide(ls_Peptide);
        double ld_BaseMass = mk_IsotopeEnvelope.massForComposition(lk_Composition);
        mk_RenderBaseMassForPeptide[ls_Peptide] = ld_BaseMass;
        
        tk_StringIntHash lk_StarAminoAcidMaxCount;
        foreach (QString ls_AminoAcid, lk_StarAminoAcidOrder)
            lk_StarAminoAcidMaxCount[ls_AminoAcid] = ls_Peptide.count(ls_AminoAcid);
        
        for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
        {
            tk_StringIntHash lk_StarAminoAcidCount;
            foreach (QString ls_AminoAcid, lk_StarAminoAcidOrder)
                lk_StarAminoAcidCount[ls_AminoAcid] = 0;
            
            while (true)
            {
                tk_IsotopeEnvelope lk_Envelope;
                if (li_Envelope == 0)
                    lk_Envelope = lightEnvelopeForPeptide(ls_Peptide);
                else
                    lk_Envelope = heavyEnvelopeForPeptide(ls_Peptide, lk_StarAminoAcidCount);
                
                tk_IsotopeEnvelope lk_EnvelopeNormalized = k_IsotopeEnvelope::normalize(lk_Envelope);
                
//                 int li_HeavyMassShift = heavyMassShiftForPeptide(ls_Peptide, lk_StarAminoAcidCount);
                
                QString ls_EnvelopeTitle = "natural";
                
                if (li_Envelope > 0)
                    ls_EnvelopeTitle = heavyEnvelopeTitle(lk_StarAminoAcidCount);
                
/*                printf("%d %9.4f %9.4f %9.4f\n", li_HeavyMassShift, lk_Envelope[li_HeavyMassShift].first, lk_Envelope[li_HeavyMassShift + 1].first, lk_Envelope[li_HeavyMassShift + 2].first);
                printf("%d %9.4f %9.4f %9.4f\n", li_HeavyMassShift, lk_Envelope[li_HeavyMassShift].second, lk_Envelope[li_HeavyMassShift + 1].second, lk_Envelope[li_HeavyMassShift + 2].second);*/
                
                for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
                {
    //                     printf("%s (E%d, %d+)\n", ls_Peptide.toStdString().c_str(), li_Envelope, li_Charge);
                    QString ls_PeptideWeight = QString("%1-%2").arg(ls_Peptide).arg(li_Envelope);
                    QString ls_PeptideWeightTitle = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Envelope).arg(ls_EnvelopeTitle);
                    QString ls_PeptideCharge = QString("%1-%2").arg(ls_Peptide).arg(li_Charge);
                    QString ls_PeptideChargeWeight = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);
                    mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].append(r_EnvelopePeaks());
                    mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].last().ms_Title = ls_EnvelopeTitle;
                    
                    if (!mk_RenderMzRangeForPeptideChargeWeight.contains(ls_PeptideChargeWeight))
                        mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight] = tk_DoublePair(1e20, -1e20);
                    
                    if (!mk_RenderIsotopeEnvelopeForPeptideWeightTitle.contains(ls_PeptideWeightTitle))
                        mk_RenderIsotopeEnvelopeForPeptideWeightTitle[ls_PeptideWeightTitle] = lk_Envelope;
                    
                    // add forbidden peak if light envelope
                    if (mb_CheckForbiddenPeak && (li_Envelope == 0))
                    {
                        double ld_Mz = (ld_BaseMass - AVERAGE_NEUTRON + md_HydrogenMass * li_Charge) / li_Charge;
                        double ld_Error = ld_Mz * (md_MassAccuracy / 1000000.0) * md_AbsenceMassAccuracyFactor;
                        double ld_MzMin = ld_Mz - ld_Error;
                        // oy, it's the forbidden peak!
                        mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].last().mk_ForbiddenIds.insert(mk_Targets.size());
                        // add dummy forbidden peak so that k_QuantifierBase::matchTargetsToPeaks 
                        // will be able to determine the mass range of the forbidden peak
                        mk_TargetMzAndIntensity[mk_Targets.size()] = tk_DoublePair(ld_Mz, 0.0);
                        mk_ForbiddenTargets.insert(mk_Targets.size());
                        mk_Targets.insert(ld_MzMin, mk_Targets.size());
    //                     /*DEBUG*/printf("%9.4f %s-%d-forbidden\n", ld_Mz, ls_Peptide.toStdString().c_str(), li_Charge);
                        mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].first = std::min<double>(ld_Mz, mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].first);
                        mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].second = std::max<double>(ld_Mz, mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].second);
                    }

                    double ld_EnvelopeBaseMass = lightMassForPeptide(ls_Peptide);
                    if (li_Envelope > 0)
                        ld_EnvelopeBaseMass += nominalMassShiftForPeptide(ls_Peptide, lk_StarAminoAcidCount);
                    double ld_EnvelopeBaseMz = (ld_EnvelopeBaseMass + md_HydrogenMass * li_Charge) / li_Charge;
                    mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].last().md_BaseMz = ld_EnvelopeBaseMz;
                    
                    // add peaks from isotope envelope
                    for (int li_PeakIndex = 0; li_PeakIndex < lk_Envelope.size(); ++li_PeakIndex)
                    {
                        bool lb_Required = false;
                        bool lb_Considered = false;
                        
                        double ld_Abundance = lk_Envelope[li_PeakIndex].first;
                        double ld_NormalizedAbundance = lk_EnvelopeNormalized[li_PeakIndex].first;
                        double ld_MassShift = lk_Envelope[li_PeakIndex].second;
                        double ld_Mz = (ld_BaseMass + ld_MassShift + md_HydrogenMass * li_Charge) / li_Charge;
                        double ld_Error = ld_Mz * (md_MassAccuracy / 1000000.0);
                        double ld_MzMin = ld_Mz - ld_Error;
                        if (mb_UseIsotopeEnvelopes)
                        {
                            lb_Required = (ld_NormalizedAbundance >= md_RequireAbundance);
                            lb_Considered = (ld_NormalizedAbundance >= md_ConsiderAbundance);
                        }
                        else
                        {
                            //int li_LightOrHeavyPeakIndex = li_PeakIndex - li_Envelope * li_HeavyMassShift;
                            //lb_Required = (li_LightOrHeavyPeakIndex >= 0) && (li_LightOrHeavyPeakIndex < mi_FixedIsotopePeakCount);
                            int li_NominalPeakIndex = (int)round(ld_BaseMass + ld_MassShift - ld_EnvelopeBaseMass);
                            lb_Required = (li_NominalPeakIndex >= 0) && (li_NominalPeakIndex < mi_FixedIsotopePeakCount) && (ld_NormalizedAbundance > 0.0001);
                        }
                        
                        if (lb_Required || lb_Considered)
                        {
                            // oy, it's a required peak!
                            if (lb_Required)
                            {
                                mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].last().mk_RequiredIds.insert(mk_Targets.size());
                                mk_PeptideChargeForRequiredTarget[mk_Targets.size()] = ls_PeptideCharge;
                            }
                            else
                                mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].last().mk_ConsideredIds.insert(mk_Targets.size());
                            
                            mk_TargetMzAndIntensity[mk_Targets.size()] = tk_DoublePair(ld_Mz, ld_Abundance);
                            if (!lb_Required)
                                mk_ConsideredTargets.insert(mk_Targets.size());
                            mk_Targets.insert(ld_MzMin, mk_Targets.size());
    //                         /*DEBUG*/printf("%9.4f %s-%d-%d-%d\n", ld_Mz, ls_Peptide.toStdString().c_str(), li_Charge, li_Envelope, li_PeakIndex);
                        }
                        
                        if (ld_NormalizedAbundance > 0.001)
                        {
                            mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].first = std::min<double>(ld_Mz, mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].first);
                            mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].second = std::max<double>(ld_Mz, mk_RenderMzRangeForPeptideChargeWeight[ls_PeptideChargeWeight].second);
                        }
                    }
                }
                // now advance the star amino acid count, and break the loop if we're finished
                if (li_Envelope == 0)
                    break;
                else
                {
                    if (lk_StarAminoAcidOrder.empty())
                        break;
                    // we're handling the heavy envelope now, increase 
                    int li_RewindUntil = -1;
                    for (int i = 0; i < lk_StarAminoAcidOrder.size(); ++i)
                    {
                        QString ls_AminoAcid = lk_StarAminoAcidOrder[i];
                        if (lk_StarAminoAcidCount[ls_AminoAcid] < lk_StarAminoAcidMaxCount[ls_AminoAcid])
                        {
                            ++lk_StarAminoAcidCount[ls_AminoAcid];
                        }
                        else
                        {
                            li_RewindUntil = i;
                            break;
                        }
                    }
                    bool lb_Finished = false;
                    if (li_RewindUntil >= 0)
                    {
                        for (int i = 0; i <= li_RewindUntil; ++i)
                            lk_StarAminoAcidCount[lk_StarAminoAcidOrder[i]] = 0;
                        if (li_RewindUntil + 1 < lk_StarAminoAcidOrder.size())
                            lk_StarAminoAcidCount[lk_StarAminoAcidOrder[li_RewindUntil + 1]] += 1;
                        else
                            lb_Finished = true;
                    }
                    if (lb_Finished)
                        break;
                }
            }
        }
    }
    
    // now check whether there are overlapping targets between individual peptide/charge combinations
    // if yes, at first remove the considered peaks
    // if this is not sufficient, remove both peptides from the list, because they are overlapping

    if (mb_CheckOverlappingPeaks)
    {
        QSet<QString> lk_OverlappingPeptideChargeItems;
        
        QSet<int> lk_IdsToBeDeleted;
        
        QMultiMap<double, int>::const_iterator lk_Iter = mk_Targets.constBegin();
        for (; lk_Iter != mk_Targets.constEnd(); ++lk_Iter)
        {
            int li_Index = lk_Iter.value();
            if (mk_ForbiddenTargets.contains(li_Index))
                continue;
            
            int li_CollidingPeaks = 0;
            
            double ld_Mz = mk_TargetMzAndIntensity[li_Index].first;
            double ld_MinMz = lk_Iter.key();
            double ld_MaxMz = 2.0 * ld_Mz - ld_MinMz;
            
            QMultiMap<double, int>::const_iterator lk_CollidingIter = lk_Iter;
            ++lk_CollidingIter;
            while (lk_CollidingIter != mk_Targets.constEnd())
            {
                double ld_OtherMz = mk_TargetMzAndIntensity[lk_CollidingIter.value()].first;
                if (ld_OtherMz > ld_MaxMz)
                    break;
                else
                    ++li_CollidingPeaks;
                ++lk_CollidingIter;
            }
            --lk_CollidingIter;
            if (li_CollidingPeaks > 0)
            {
                // there are colliding peaks, remove all candidate peaks 
                // and see if only one required peak remains
                ++lk_CollidingIter;
                QSet<int> lk_RequiredTargets;
                QMultiMap<double, int>::const_iterator lk_Temp = lk_Iter;
                for (; lk_Temp != lk_CollidingIter; ++lk_Temp)
                {
                    int li_Index = lk_Temp.value();
                    if (mk_ForbiddenTargets.contains(li_Index))
                        continue;
                    
                    if (mk_ConsideredTargets.contains(li_Index))
                    {
                        // it's a considered target, remove it
                        lk_IdsToBeDeleted << li_Index;
                    }
                    else
                        lk_RequiredTargets << li_Index;
                }
                if (lk_RequiredTargets.size() > 1)
                {
                    // we still have more than one required peak here,
                    // mark all affected required peaks for deletetion,
                    // which will lead to non-handling of the appropriate peptides
                    // because at least one required peak could not be found,
                    // print a message to the screen later
                    foreach (int li_Id, lk_RequiredTargets)
                    {
                        lk_OverlappingPeptideChargeItems << mk_PeptideChargeForRequiredTarget[li_Id];
                        lk_IdsToBeDeleted << li_Index;
                    }
                }
            }
        }
        
        // delete items to be deleted by re-creating mk_Targets
        QMultiMap<double, int> lk_NewTargets;
        lk_Iter = mk_Targets.constBegin();
        for (; lk_Iter != mk_Targets.constEnd(); ++lk_Iter)
        {
            int li_Index = lk_Iter.value();
            if (!lk_IdsToBeDeleted.contains(li_Index))
                lk_NewTargets.insert(lk_Iter.key(), lk_Iter.value());
        }
        mk_Targets = lk_NewTargets;
        foreach (int li_Id, lk_IdsToBeDeleted)
            mk_TargetMzAndIntensity.remove(li_Id);
        
        if (!lk_OverlappingPeptideChargeItems.empty())
        {
            QHash<QString, int> lk_PeptideRemoveCount;
            foreach (QString ls_PeptideCharge, lk_OverlappingPeptideChargeItems)
            {
                QString ls_Peptide = ls_PeptideCharge.split("-").first();
                if (!lk_PeptideRemoveCount.contains(ls_Peptide))
                    lk_PeptideRemoveCount[ls_Peptide] = 0;
                ++lk_PeptideRemoveCount[ls_Peptide];
            }
            int li_RemovedAllCount = 0;
            int li_RemovedSomeCount = 0;
            foreach (QString ls_Peptide, lk_PeptideRemoveCount.keys())
            {
                if (lk_PeptideRemoveCount[ls_Peptide] == (mi_MaxCharge - mi_MinCharge + 1))
                    ++li_RemovedAllCount;
                else
                    ++li_RemovedSomeCount;
            }
            QStringList lk_Message;
            if (li_RemovedAllCount > 0)
                lk_Message << QString("completely ignoring %1 peptide%2").arg(li_RemovedAllCount).arg(li_RemovedAllCount != 1 ? "s" : "");
            if (li_RemovedSomeCount > 0)
                lk_Message << QString("ignoring some charge states of %1 peptide%2").arg(li_RemovedSomeCount).arg(li_RemovedSomeCount != 1 ? "s" : "");
            QString ls_Message = lk_Message.join(" and ");
            printf("Warning: %s because of overlapping required peaks.\n", ls_Message.toStdString().c_str());
        }
    }
    
    
    printf("done.\n");
    
    fflush(stdout);

/*    for (int i = 0; i < mk_Targets.size(); ++i)
        printf("%3d %9.4f %3d\n", i, mk_Targets.keys()[i], mk_Targets.values()[i]);*/
    
    if (mk_pCsvStream)
    {
        *(mk_pCsvStream.data()) << "Filename,Scan id,Peptide,Charge,Amount light,Amount heavy,Retention time";
        if (mb_UseIsotopeEnvelopes)
            *(mk_pCsvStream.data()) << ",Mean error light,Max error light,Mean error heavy,Max error heavy";
        *(mk_pCsvStream.data()) << "\n";
    }
    
    
    if (mk_pXhtmlStream)
    {
        QFile lk_File(":res/qtrace-xhtml-header.xhtml.part");
        lk_File.open(QIODevice::ReadOnly);
        QByteArray lk_Content = lk_File.readAll();
        *(mk_pXhtmlStream.data()) << QString(lk_Content);
        *(mk_pXhtmlStream.data()) << "<thead><tr><th>Filename</th><th>Scan</th><th>Peptide</th><th>Charge</th><th>Amount light</th><th>Amount heavy</th><th>Retention time</th>";
/*        if (mb_UseIsotopeEnvelopes)
            *(mk_pXhtmlStream.data()) << "<th>Mean error light</th><th>Max error light</th><th>Mean error heavy</th><th>Max error heavy</th>";*/
        *(mk_pXhtmlStream.data()) << "</tr></thead>\n<tbody>\n";

        lk_File.close();
    }
    
    // parse all spectra files
    foreach (QString ls_Path, mk_SpectraFiles)
    {
        ms_CurrentSpectraFile = QFileInfo(ls_Path).baseName();

        this->parseFile(ls_Path);
        
        if (!mb_Quiet)
            printf(" done.\n");
    }
    
    if (mk_pXhtmlStream)
    {
        *(mk_pXhtmlStream.data()) << "</tbody>\n";
        QFile lk_File(":res/qtrace-xhtml-footer.xhtml.part");
        lk_File.open(QIODevice::ReadOnly);
        QByteArray lk_Content = lk_File.readAll();
        *(mk_pXhtmlStream.data()) << QString(lk_Content);
        lk_File.close();
    }
}


void k_Quantifier::handleScan(r_Scan& ar_Scan, bool& ab_Continue)
{
    ab_Continue = true;
    
    if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
    {
        printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
        return;
    }

    // find all peaks in this spectrum
    QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum, md_MinSnr);
    QList<double> lk_AllPeakMz;
    foreach (r_Peak lr_Peak, lk_AllPeaks)
        lk_AllPeakMz.append(lr_Peak.md_PeakMz);
    
    QHash<int, int> lk_PreMatches = matchTargetsToPeaks(lk_AllPeakMz, mk_Targets, mk_TargetMzAndIntensity, mk_ForbiddenTargets);
    QHash<int, int> lk_Matches;
    foreach (int li_Key, lk_PreMatches.keys())
        lk_Matches[mk_Targets.values()[li_Key]] = lk_PreMatches[li_Key];
    
    QSet<int> lk_MatchedTargetIds = lk_Matches.keys().toSet();
    
    foreach (QString ls_Peptide, mk_Peptides)
    {
        for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
        {
            r_ScanQuantitationResult lr_ScanResult;

            // try to quantify (ls_Peptide) at (li_Charge) in the current scan,
            // and quantify light and heavy peptides independently,
            // except if a forbidden peak should be checked for and is actually there 
            
            // if the forbidden peak is there, don't use this at all!
            bool lb_SeenForbiddenPeak = false;
            
            QHash<QString, double> lk_Amounts;
            QHash<QString, double> lk_FitError;
            QHash<QString, double> lk_IndividualFitError;
            QHash<QString, bool> lk_Good;
            
            bool lb_GoodEnvelope_[2] = {true, true};
            double ld_AmountsEnvelope_[2] = {0.0, 0.0};
            QList<double> lk_FitErrors_[2] = {QList<double>(), QList<double>()};
            
            for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
            {
                QString ls_PeptideChargeWeightKey = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);
                
                foreach (r_EnvelopePeaks lr_Peaks, mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey])
                {
                    QString ls_EnvelopeTitle = lr_Peaks.ms_Title;
                    QString ls_WeightTitle = QString("%1-%2").arg(li_Envelope).arg(ls_EnvelopeTitle);
                    
                    lk_Amounts[ls_WeightTitle] = 0.0;
                    lk_FitError[ls_WeightTitle] = 0.0;
                    lk_IndividualFitError[ls_WeightTitle] = 0.0;
                    lk_Good[ls_WeightTitle] = true;
                    
                    // don't quantify this at all if the forbidden peak was there
                    if ((li_Envelope == 0) && (mb_CheckForbiddenPeak && (!(lr_Peaks.mk_ForbiddenIds & lk_MatchedTargetIds).empty())))
                    {
                        lb_SeenForbiddenPeak = true;
                        break;
                    }
                    
                    // skip this if not all required light and heavy peaks are there
                    if ((lk_MatchedTargetIds & lr_Peaks.mk_RequiredIds).size() != lr_Peaks.mk_RequiredIds.size())
                    {
                        lk_Good[ls_WeightTitle] = false;
                    }
                    else
                    {
                        if (mb_UseIsotopeEnvelopes)
                        {
                            // extreme isotope envelope least squares fitting action...
                            QSet<int> lk_AllTargets = 
                                lr_Peaks.mk_RequiredIds | 
                                lr_Peaks.mk_ConsideredIds;
                                
                            QList<tk_DoublePair> lk_Pairs = QList<tk_DoublePair>();
                            
                            foreach (int li_Id, lk_AllTargets)
                            {
                                if (lk_MatchedTargetIds.contains(li_Id))
                                {
                                    double ld_TargetIntensity = mk_TargetMzAndIntensity[li_Id].second;
                                    double ld_PeakIntensity = lk_AllPeaks[lk_Matches[li_Id]].md_PeakIntensity;
                                    lk_Pairs << tk_DoublePair(ld_TargetIntensity, ld_PeakIntensity);
                                }
                            }
//                             printf("pairs: %d / %s %d\n", lk_Pairs.size(), ls_Peptide.toStdString().c_str(), li_Charge);
                            
                            double ld_FitFactor = 0.0;
                            if (lk_Pairs.size() > 0)
                            {
                                QList<double> lk_ThisFitErrors;
                                leastSquaresFit(lk_Pairs, &ld_FitFactor, &lk_ThisFitErrors);
                                double ld_FitError = 0.0;
                                double ld_IndividualFitError = 0.0;
                                lk_FitErrors_[li_Envelope] += lk_ThisFitErrors;
                                foreach (double ld_Error, lk_ThisFitErrors)
                                    ld_FitError += ld_Error;
                                ld_FitError /= lk_ThisFitErrors.size();
                                
                                if (ld_FitError > md_MaxFitError)
                                    lk_Good[ls_WeightTitle] = false;
                                lk_FitError[ls_WeightTitle] = ld_FitError;
                                lr_ScanResult.mk_ProfileScale[ls_WeightTitle] = ld_FitFactor;
                                lk_IndividualFitError[ls_WeightTitle] = ld_IndividualFitError;
/*                                foreach (tk_DoublePair lk_Pair, lk_Pairs)
                                    printf("%9.4f - %9.4f\n", lk_Pair.first * ld_FitFactor, lk_Pair.second);*/
                            }
                            else
                                lk_Good[ls_WeightTitle] = false;
                            
                            if (lk_Good[ls_WeightTitle])
                                lk_Amounts[ls_WeightTitle] += ld_FitFactor;
                            lr_ScanResult.mk_Good[ls_WeightTitle] = lk_Good[ls_WeightTitle];
                        }
                        else
                        {
                            // just add the peak heights, that's all
                            foreach (int li_Id, lr_Peaks.mk_RequiredIds)
                                lk_Amounts[ls_WeightTitle] += lk_AllPeaks[lk_Matches[li_Id]].md_PeakIntensity;
                        }
                    }
                    ld_AmountsEnvelope_[li_Envelope] += lk_Amounts[ls_WeightTitle];
                    if (!lk_Good[ls_WeightTitle])
                        lb_GoodEnvelope_[li_Envelope] = false;
                }
                
                if (lb_SeenForbiddenPeak)
                    break;
            }
            
            // set envelope amounts to 0.0 if not good
            for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
            {
                if (!lb_GoodEnvelope_[li_Envelope])
                    ld_AmountsEnvelope_[li_Envelope] = 0.0;
            }
            
            if (!lb_SeenForbiddenPeak)
            {
                if (lb_GoodEnvelope_[0] || lb_GoodEnvelope_[1])
                {
                    if (mk_pCsvStream)
                    {
                        *(mk_pCsvStream.data()) << QString("\"%1\",\"%2\",%3,%4,%5,%6,%7")
                            .arg(ms_CurrentSpectraFile)
                            .arg(ar_Scan.ms_Id)
                            .arg(ls_Peptide)
                            .arg(li_Charge)
                            .arg(ld_AmountsEnvelope_[0], 1, 'f', 4)
                            .arg(ld_AmountsEnvelope_[1], 1, 'f', 4)
                            .arg(ar_Scan.md_RetentionTime, 1, 'f', 2);
                        if (mb_UseIsotopeEnvelopes)
                        {
                            for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
                            {
                                double ld_MeanError = 0.0;
                                double ld_MaxError = 0.0;
                                foreach (double ld_Error, lk_FitErrors_[li_Envelope])
                                {
                                    ld_MeanError += ld_Error;
                                    ld_MaxError = std::max<double>(ld_MaxError, ld_Error);
                                }
                                ld_MeanError /= lk_FitErrors_[li_Envelope].size();
                                *(mk_pCsvStream.data()) << QString(",%1,%2") 
                                    .arg(ld_MeanError, 1, 'f', 3)
                                    .arg(ld_MaxError, 1, 'f', 3);
                            }
                        }
                        *(mk_pCsvStream.data()) << "\n";
                    }
                    if (mk_pXhtmlStream)
                    {
                        *(mk_pXhtmlStream.data()) << QString("\n<!-- BEGIN PEPTIDE %1 CHARGE %2 BAND %3 SCAN %4 -->\n").arg(ls_Peptide).arg(li_Charge).arg(ms_CurrentSpectraFile).arg(ar_Scan.ms_Id);
                        *(mk_pXhtmlStream.data()) << QString("<tr><td>%1</td><td>%2</td><td>%3</td><td>%4</td><td>%5</td><td>%6</td><td>%7</td>")
                            .arg(ms_CurrentSpectraFile)
                            .arg(ar_Scan.ms_Id)
                            .arg(ls_Peptide)
                            .arg(li_Charge)
                            .arg(ld_AmountsEnvelope_[0], 1, 'f', 4)
                            .arg(ld_AmountsEnvelope_[1], 1, 'f', 4)
                            .arg(ar_Scan.md_RetentionTime, 1, 'f', 2);
                            
/*                        if (mb_UseIsotopeEnvelopes)
                        {
                            for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
                            {
                                double ld_MeanError = 0.0;
                                double ld_MaxError = 0.0;
                                foreach (double ld_Error, lk_FitErrors_[li_Envelope])
                                {
                                    ld_MeanError += ld_Error;
                                    ld_MaxError = std::max<double>(ld_MaxError, ld_Error);
                                }
                                ld_MeanError /= lk_FitErrors_[li_Envelope].size();
                                *(mk_pXhtmlStream.data()) << QString("<td>%1%</td><td>%2%</td>")
                                    .arg(ld_MeanError * 100.0, 1, 'f', 1)
                                    .arg(ld_MaxError * 100.0, 1, 'f', 1);
                            }
                        }*/
/*                        if (mb_UseIsotopeEnvelopes)
                            *(mk_pXhtmlStream.data()) << QString("<td>%1%</td><td>%2%</td><td>%3%</td><td>%4%</td>")
                                .arg(ld_FitError_[0] * 100.0, 1, 'f', 1)
                                .arg(ld_IndividualFitError_[0] * 100.0, 1, 'f', 1)
                                .arg(ld_FitError_[1] * 100.0, 1, 'f', 1)
                                .arg(ld_IndividualFitError_[1] * 100.0, 1, 'f', 1);*/
                            
                        *(mk_pXhtmlStream.data()) << "</tr>\n";
                            
                        lr_ScanResult.ms_Peptide = ls_Peptide;
                        lr_ScanResult.mi_Charge = li_Charge;
                        lr_ScanResult.md_AmountUnlabeled = ld_AmountsEnvelope_[0];
                        lr_ScanResult.md_AmountLabeled = ld_AmountsEnvelope_[1];
                        
                        lr_ScanResult.md_MinMz = mk_RenderMzRangeForPeptideChargeWeight[QString("%1-%2-0").arg(ls_Peptide).arg(li_Charge)].first;
                        lr_ScanResult.md_MaxMz = mk_RenderMzRangeForPeptideChargeWeight[QString("%1-%2-1").arg(ls_Peptide).arg(li_Charge)].second;
                            
                        QString ls_Svg = this->renderScanAsSvg(ar_Scan, lr_ScanResult, lk_AllPeaks, lk_Matches);
                        ls_Svg.remove(QRegExp("<\\?xml.+\\?>"));
                        ls_Svg.replace(QRegExp("width=\\\"[^\\\"]*\\\"\\s+height=\\\"[^\\\"]*\\\""), "width='950' height='238'");
                        *(mk_pXhtmlStream.data()) << "<div style='background-color: #fff;' width='950' height='238'>";
                        *(mk_pXhtmlStream.data()) << ls_Svg;
                        *(mk_pXhtmlStream.data()) << "</div>" << endl;
                        *(mk_pXhtmlStream.data()) << QString("\n<!-- END PEPTIDE %1 CHARGE %2 BAND %3 SCAN %4 -->\n").arg(ls_Peptide).arg(li_Charge).arg(ms_CurrentSpectraFile).arg(ar_Scan.ms_Id);
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
            else
            {
                printf("Error: Invalid argument specified: %s.\n", ls_Key.toStdString().c_str());
                exit(1);
            }
        }
    }
    if (!ak_Arguments.empty())
    {
        printf("Error: Invalid argument specified: %s.\n", ak_Arguments.first().toStdString().c_str());
        exit(1);
    }
}
