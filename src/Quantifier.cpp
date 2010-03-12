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
                    QString ls_PeptideChargeWeight = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);
                    mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight] = r_EnvelopePeaks();
                    
                    // add forbidden peak if light envelope
                    if (mb_CheckForbiddenPeak && (li_Envelope == 0))
                    {
                        double ld_Mz = (ld_BaseMass -md_HydrogenMass + md_HydrogenMass * li_Charge) / li_Charge;
                        // oy, it's the forbidden peak!
                        mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].mk_ForbiddenIds.insert(mk_Targets.size());
                        mk_Targets.insert(ld_Mz, mk_Targets.size());
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
                        }
                        else if (ld_NormalizedAbundance >= md_ConsiderAbundance)
                        {
                            // oy, it's a peak to be considered!
                            mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeight].mk_ConsideredIds.insert(mk_Targets.size());
                            mk_TargetIntensity[mk_Targets.size()] = ld_Abundance;
                            mk_Targets.insert(ld_Mz, mk_Targets.size());
                        }
                    }
                }
            }
        }
    }
    
    if (mk_pCsvStream)
        *(mk_pCsvStream.get_Pointer()) << "filename,scan id,peptide,amount light,amount heavy,retention time,charge\n";
    
    // parse all spectra files
    foreach (QString ls_Path, mk_SpectraFiles)
    {
        ms_CurrentSpectraFile = QFileInfo(ls_Path).baseName();

        this->parseFile(ls_Path);
        
        if (!mb_Quiet)
            printf(" done.\n");
    }
    /*
	mk_Peptides = ak_Peptides;
	QSet<QString> lk_PeptidesActuallySearchedFor;
    
    typedef QPair<double, QString> tk_DoubleStringPair;
    QList<tk_DoubleStringPair> lk_TempList;
    // determine all peptide isotope envelopes and target m/z values
    foreach (QString ls_Peptide, mk_Peptides)
    {
        // check whether we have to skip this peptide
        // :TODO: skip a peptide if no mass shift!!!!
        lk_PeptidesActuallySearchedFor << ls_Peptide;
        
        QHash<QString, int> lk_Composition = compositionForPeptide(ls_Peptide);
        tk_IsotopeEnvelope lk_UnchargedIsotopeEnvelope = mk_IsotopeEnvelope.isotopeEnvelopeForComposition(lk_Composition);
        double ld_UnchargedPeptideMass = mk_IsotopeEnvelope.massForComposition(lk_Composition);
        tk_IsotopeEnvelope lk_IsotopeEnvelope = lk_UnchargedIsotopeEnvelope;
        tk_IsotopeEnvelope lk_UnchargedIsotopeEnvelopeHeavy = heavyEnvelopeForPeptide(ls_Peptide);
        tk_IsotopeEnvelope lk_IsotopeEnvelopeHeavy = lk_UnchargedIsotopeEnvelopeHeavy;
        for (int li_Charge = 1; li_Charge <= mi_MaxCharge; ++li_Charge)
        {
            QString ls_PeptideChargeKey = QString("%1-%2").arg(ls_Peptide).arg(li_Charge);
            double ld_Mass = ld_UnchargedPeptideMass + mk_IsotopeEnvelope.mk_BaseIsotopeMass["H"] * li_Charge;
            lk_IsotopeEnvelope = mk_IsotopeEnvelope.add(lk_IsotopeEnvelope, mk_IsotopeEnvelope.mk_ElementEnvelopes["H"].first());
            lk_IsotopeEnvelopeHeavy = mk_IsotopeEnvelope.add(lk_IsotopeEnvelopeHeavy, mk_IsotopeEnvelope.mk_ElementEnvelopes["H"].first());
            
            if (li_Charge >= mi_MinCharge && li_Charge <= mi_MaxCharge)
            {
//                  fprintf(stderr, "%s (%d+)\n", ls_Peptide.toStdString().c_str(), li_Charge);
                mk_UnlabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey] = lk_IsotopeEnvelope;
                tk_IsotopeEnvelope lk_NormalizedEnvelope = mk_IsotopeEnvelope.normalize(lk_IsotopeEnvelope);

                // store unlabeled isotope envelope
                bool lb_SeenRequired = false;
                for (int li_Isotope = 0; li_Isotope < lk_NormalizedEnvelope.size(); ++li_Isotope)
                {
                    if (lk_NormalizedEnvelope[li_Isotope].first >= md_ConsiderAbundance)
                    {
                        double ld_Mz = (ld_Mass + lk_NormalizedEnvelope[li_Isotope].second) / li_Charge;
                        // peptide-charge-label-isotope
                        QString ls_Key = QString("%1-%2-0-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Isotope);
                        lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));
                        
                        if (lk_NormalizedEnvelope[li_Isotope].first >= md_RequireAbundance)
                        {
                            mk_UnlabeledRequiredTargetMzForPeptideCharge[ls_PeptideChargeKey] << ls_Key;
                            lb_SeenRequired = true;
                        }
                        else
                        {
                            if (!lb_SeenRequired)
                                mk_UnlabeledConsideredLeftTargetMzForPeptideCharge[ls_PeptideChargeKey].insert(0, ls_Key);
                            else
                                mk_UnlabeledConsideredRightTargetMzForPeptideCharge[ls_PeptideChargeKey].append(ls_Key);
                        }

//                        fprintf(stderr, "A+%d %c %9.5f / %9.5f\n", li_Isotope, 
//                                lk_NormalizedEnvelope[li_Isotope].first >= REQUIRE_ABUNDANCE ? '*' : lk_NormalizedEnvelope[li_Isotope].first >= CONSIDER_ABUNDANCE ? '+' : ' ', 
//                                lk_NormalizedEnvelope[li_Isotope].first, ld_Mz);
                    }
                }
                
                // store labeled isotope envelope
                mk_LabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey] = lk_IsotopeEnvelopeHeavy;
                tk_IsotopeEnvelope lk_NormalizedEnvelopeHeavy = mk_IsotopeEnvelope.normalize(lk_IsotopeEnvelopeHeavy);
                lb_SeenRequired = false;
                for (int li_Isotope = 0; li_Isotope < lk_NormalizedEnvelopeHeavy.size(); ++li_Isotope)
                {
                    if (lk_NormalizedEnvelopeHeavy[li_Isotope].first >= md_ConsiderAbundance)
                    {
                        double ld_Mz = (ld_Mass + lk_NormalizedEnvelopeHeavy[li_Isotope].second) / li_Charge;
                        // peptide-charge-label-isotope
                        QString ls_Key = QString("%1-%2-1-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Isotope);
                        lk_TempList.push_back(tk_DoubleStringPair(ld_Mz, ls_Key));

                        if (lk_NormalizedEnvelopeHeavy[li_Isotope].first >= md_RequireAbundance)
                        {
                            mk_LabeledRequiredTargetMzForPeptideCharge[ls_PeptideChargeKey] << ls_Key;
                            lb_SeenRequired = true;
                        }
                        else
                        {
                            if (!lb_SeenRequired)
                                mk_LabeledConsideredLeftTargetMzForPeptideCharge[ls_PeptideChargeKey].insert(0, ls_Key);
                            else
                                mk_LabeledConsideredRightTargetMzForPeptideCharge[ls_PeptideChargeKey].append(ls_Key);
                        }

//                        fprintf(stderr, "A*+%d %c %9.5f / %9.5f\n", li_Isotope, 
//                                lk_NormalizedEnvelopeN15[li_Isotope].first > REQUIRE_ABUNDANCE ? '*' : lk_NormalizedEnvelopeN15[li_Isotope].first > CONSIDER_ABUNDANCE ? '+' : ' ', 
//                                lk_NormalizedEnvelopeN15[li_Isotope].first, ld_Mz);
                    }
                }
            }
        }
    }
    
	if (!mb_Quiet)
	{
		if (lk_PeptidesActuallySearchedFor.empty())
			printf("No appropriate peptides left for the search, skipping input files...");
        else
        {
            printf("Searching for %d peptide%s in %d file%s, trying charge states %d to %d, requiring a SNR of %1.2f.\n",
                lk_PeptidesActuallySearchedFor.size(), ak_Peptides.size() != 1 ? "s" : "",
                ak_SpectraFiles.size(), ak_SpectraFiles.size() != 1 ? "s" : "",
                mi_MinCharge, mi_MaxCharge, md_MinSnr);
            printf("Requiring peaks down to %1.3f, considering down to %1.3f, allowing for max fit error of %1.3f.\n", md_RequireAbundance, md_ConsiderAbundance, md_MaxFitError);
            QStringList lk_ScanTypes;
            if (me_ScanType & r_ScanType::MS1)
                lk_ScanTypes << "Full";
            if (me_ScanType & r_ScanType::SIM)
                lk_ScanTypes << "SIM";

            printf("Looking in %s scans, mass accuracy is %1.2f ppm, %s the forbidden peak.\n",
                lk_ScanTypes.join("/").toStdString().c_str(), md_MassAccuracy, mb_CheckForbiddenPeak ? "checking for" : "ignoring");
        }
	}
	
    // :TODO: use a map here
	//qSort(lk_TempList.begin(), lk_TempList.end(), sortByMz);

	mk_AllTargetMasses = QList<double>();
	mk_TargetMzIndex = QHash<QString, int>();

	foreach (tk_DoubleStringPair lk_Pair, lk_TempList)
	{
		mk_TargetMzIndex[lk_Pair.second] = mk_AllTargetMasses.size();
		mk_AllTargetMasses.push_back(lk_Pair.first);
	}
	
	if (mk_CsvOutStream.device())
		mk_CsvOutStream << "Filename,Scan id,Peptide,Amount light,Amount heavy,Retention time,Charge,Filter line,SNR" << endl;
	
	if (mk_XhtmlOutStream.device())
	{
		QFile lk_File(":res/qtrace-xhtml-header.xhtml.part");
		lk_File.open(QIODevice::ReadOnly);
		QByteArray lk_Content = lk_File.readAll();
		mk_XhtmlOutStream << QString(lk_Content);
        mk_XhtmlOutStream << QString("<table>\n<tr><th>Filename</th><th>Scan id</th><th>Peptide</th><th>Amount light</th><th>Amount heavy</th><th>Retention time</th><th>Charge</th><th>Filter line</th><th>SNR</th></tr>\n");
		lk_File.close();
	}
		
	// parse all bands
	foreach (QString ls_Path, ak_SpectraFiles)
	{
		ms_CurrentSpot = QFileInfo(ls_Path).baseName();
		
		// parse spot
		if (!lk_PeptidesActuallySearchedFor.empty())
			this->parseFile(ls_Path);
		
		if (!mb_Quiet)
			printf(" done.\n");
	}

	if (mk_XhtmlOutStream.device())
	{
        mk_XhtmlOutStream << QString("</table>\n");
		QFile lk_File(":res/qtrace-xhtml-footer.xhtml.part");
		lk_File.open(QIODevice::ReadOnly);
		QByteArray lk_Content = lk_File.readAll();
		mk_XhtmlOutStream << QString(lk_Content);
		lk_File.close();
	}
    */
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
        lk_AllPeakMz.append(lr_Peak.md_PeakMz);
    
    QHash<int, int> lk_Matches = matchTargetsToPeaks(lk_AllPeakMz, mk_Targets, md_MassAccuracy);
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
                bool lb_Good_[2] = {true, true};
                
                for (int li_Envelope = 0; li_Envelope < 2; ++li_Envelope)
                {
                    QString ls_PeptideChargeWeightKey = QString("%1-%2-%3").arg(ls_Peptide).arg(li_Charge).arg(li_Envelope);

                    // don't quantify this at all if the forbidden peak was there
                    if (li_Envelope == 0 && (mb_CheckForbiddenPeak && (!(mk_TargetsForPeptideChargeWeight[ls_PeptideChargeWeightKey].mk_ForbiddenIds & lk_MatchedTargetIds).empty())))
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
                            lb_Good_[0] = lb_Good_[1] = false;
                            break;
                        }
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
                            .arg(ld_FitFactor_[0], 1, 'f', 4)
                            .arg(ld_FitFactor_[1], 1, 'f', 4)
                            .arg(ar_Scan.md_RetentionTime)
                            .arg(li_Charge);
                    }
                }
            }
        }
    }
    
    /*

    QList<double> lk_AllPeaksMz;
    foreach (r_Peak lr_Peak, lk_AllPeaks)
        lk_AllPeaksMz << lr_Peak.md_PeakMz;
    
    QHash<int, int> lk_PeakForTargetMz = matchTargetsToPeaks(lk_AllPeaksMz, mk_AllTargetMasses, md_MassAccuracy);
	
//	foreach (QString s, mk_TargetMzIndex.keys())
// 	{
// 		double ld_ScanMz = lk_AllPeaks[lk_PeakForTargetMz[mk_TargetMzIndex[s]]].md_PeakMz;
// 		double ld_TargetMz = mk_AllTargetMasses[mk_TargetMzIndex[s]];
// 		double ld_Error = fabs(ld_ScanMz - ld_TargetMz) / ld_TargetMz * 1000000.0;
// 		double ld_Snr = lk_AllPeaks[lk_PeakForTargetMz[mk_TargetMzIndex[s]]].md_Snr;
// 		printf("%s: %1.6f - %1.6f (%1.2f ppm, SNR %1.2f)\n", s.toStdString().c_str(), ld_ScanMz, ld_TargetMz, ld_Error, ld_Snr);
// 	}
	
	// discard all target m/z matches that are bogus because
	// they are not within the specified mass accuracy
	QHash<int, int> lk_IncludePeakForTargetMz;
	QHash<int, int> lk_ExcludePeakForTargetMz;
	foreach (int li_TargetMzIndex, lk_PeakForTargetMz.keys())
	{
		double ld_TargetMz = mk_AllTargetMasses[li_TargetMzIndex];
		r_Peak lr_Peak = lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];
		double ld_PeakMz = lr_Peak.md_PeakMz;
		double ld_MzError = fabs(ld_TargetMz - ld_PeakMz);
		double ld_MaxMzError = ld_TargetMz * md_MassAccuracy / 1000000.0;
		if (ld_MzError <= ld_MaxMzError)
			lk_IncludePeakForTargetMz[li_TargetMzIndex] = lk_PeakForTargetMz[li_TargetMzIndex];
		ld_MaxMzError = ld_TargetMz * (md_MassAccuracy * ABSENCE_MASS_ACCURACY_FACTOR) / 1000000.0;
		if (ld_MzError <= ld_MaxMzError)
			lk_ExcludePeakForTargetMz[li_TargetMzIndex] = lk_PeakForTargetMz[li_TargetMzIndex];
	}

    QSet<QString> lk_AvailableIncludePeaks;
    foreach (QString ls_Key, mk_TargetMzIndex.keys())
        if (lk_IncludePeakForTargetMz.contains(mk_TargetMzIndex[ls_Key]))
            lk_AvailableIncludePeaks << ls_Key;
    
    foreach (QString ls_Peptide, mk_Peptides)
    {
        for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
        {
            QString ls_PeptideChargeKey = QString("%1-%2").arg(ls_Peptide).arg(li_Charge);
            QSet<QString> lk_InterestingKeys;
            lk_InterestingKeys |= mk_UnlabeledRequiredTargetMzForPeptideCharge[ls_PeptideChargeKey];
            lk_InterestingKeys |= mk_UnlabeledConsideredLeftTargetMzForPeptideCharge[ls_PeptideChargeKey].toSet();
            lk_InterestingKeys |= mk_UnlabeledConsideredRightTargetMzForPeptideCharge[ls_PeptideChargeKey].toSet();
            lk_InterestingKeys |= mk_LabeledRequiredTargetMzForPeptideCharge[ls_PeptideChargeKey];
            lk_InterestingKeys |= mk_LabeledConsideredLeftTargetMzForPeptideCharge[ls_PeptideChargeKey].toSet();
            lk_InterestingKeys |= mk_LabeledConsideredRightTargetMzForPeptideCharge[ls_PeptideChargeKey].toSet();
            double ld_MinMz = mk_AllTargetMasses[mk_TargetMzIndex[lk_InterestingKeys.toList().first()]];
            double ld_MaxMz = ld_MinMz;
            foreach (QString ls_Key, lk_InterestingKeys)
            {
                double ld_Mz = mk_AllTargetMasses[mk_TargetMzIndex[ls_Key]];
                ld_MinMz = std::min<double>(ld_MinMz, ld_Mz);
                ld_MaxMz = std::max<double>(ld_MaxMz, ld_Mz);
            }
            double ld_MzBorder = (ld_MaxMz - ld_MinMz) * 0.01;
            ld_MinMz -= 1.0 / li_Charge;
            ld_MinMz -= ld_MzBorder;
            ld_MaxMz += ld_MzBorder;
            QSet<QString> lk_RequiredUnlabeledKeys = mk_UnlabeledRequiredTargetMzForPeptideCharge[ls_PeptideChargeKey];
            QSet<QString> lk_RequiredLabeledKeys = mk_LabeledRequiredTargetMzForPeptideCharge[ls_PeptideChargeKey];
            bool lb_FoundLightEnvelope = (lk_RequiredUnlabeledKeys & lk_AvailableIncludePeaks).size() == lk_RequiredUnlabeledKeys.size();
            bool lb_FoundHeavyEnvelope = (lk_RequiredLabeledKeys & lk_AvailableIncludePeaks).size() == lk_RequiredLabeledKeys.size();
            if (lb_FoundLightEnvelope || lb_FoundHeavyEnvelope)
            {
                // oy, we found all required peaks of at least one state!
//                 printf("%s %d+ %d %d\n", ls_Peptide.toStdString().c_str(),
//                        li_Charge, lb_FoundLightEnvelope, lb_FoundHeavyEnvelope);
                
                QSet<QString> lk_UnlabeledKeys;
                QSet<QString> lk_LabeledKeys;
                
                // if we found the light envelope, expand to left and right
                if (lb_FoundLightEnvelope)
                {
                    lk_UnlabeledKeys = lk_RequiredUnlabeledKeys;
                    foreach (QString ls_Key, mk_UnlabeledConsideredLeftTargetMzForPeptideCharge[ls_PeptideChargeKey])
                    {
                        if (lk_AvailableIncludePeaks.contains(ls_Key))
                            lk_UnlabeledKeys << ls_Key;
                        else
                            break;
                    }
                    foreach (QString ls_Key, mk_UnlabeledConsideredRightTargetMzForPeptideCharge[ls_PeptideChargeKey])
                    {
                        if (lk_AvailableIncludePeaks.contains(ls_Key))
                            lk_UnlabeledKeys << ls_Key;
                        else
                            break;
                    }
                }
                
                // if we found the heavy envelope, expand to left and right
                if (lb_FoundHeavyEnvelope)
                {
                    lk_LabeledKeys = lk_RequiredLabeledKeys;
                    foreach (QString ls_Key, mk_LabeledConsideredLeftTargetMzForPeptideCharge[ls_PeptideChargeKey])
                    {
                        if (lk_AvailableIncludePeaks.contains(ls_Key))
                            lk_LabeledKeys << ls_Key;
                        else
                            break;
                    }
                    foreach (QString ls_Key, mk_LabeledConsideredRightTargetMzForPeptideCharge[ls_PeptideChargeKey])
                    {
                        if (lk_AvailableIncludePeaks.contains(ls_Key))
                            lk_LabeledKeys << ls_Key;
                        else
                            break;
                    }
                }
                
                r_ScanQuantitationResult lr_ScanResult;
                lr_ScanResult.mb_IsGood = true;
                lr_ScanResult.ms_Peptide = ls_Peptide;
                lr_ScanResult.mi_Charge = li_Charge;
                lr_ScanResult.md_Snr = 10000.0;
                lr_ScanResult.md_RetentionTime = ar_Scan.md_RetentionTime;
                lr_ScanResult.md_MinMz = ld_MinMz;
                lr_ScanResult.md_MaxMz = ld_MaxMz;
                lr_ScanResult.md_AmountUnlabeled = 0.0;
                lr_ScanResult.md_AmountLabeled = 0.0;
                lr_ScanResult.md_UnlabeledError = 0.0;
                lr_ScanResult.md_LabeledError = 0.0;
                
                if (!lk_UnlabeledKeys.empty())
                {
                    QList<tk_DoublePair> lk_MatchValues;
                    foreach (QString ls_Key, lk_UnlabeledKeys)
                    {
                        int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
                        r_Peak lr_Peak = lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];
//                         if (me_AmountEstimation == r_AmountEstimation::Area)
//                             lr_ScanResult.md_AmountUnlabeled += lr_Peak.md_PeakArea;
//                         else if (me_AmountEstimation == r_AmountEstimation::Intensity)
                        // :TODO: now this always uses peak height, make this an option (area/intensity)
                        lr_ScanResult.md_AmountUnlabeled += lr_Peak.md_PeakIntensity;
                        lr_ScanResult.md_Snr = std::min<double>(lr_Peak.md_Snr, lr_ScanResult.md_Snr);
                        lr_ScanResult.mk_UnlabeledPeaks.append(lr_Peak);
                        int li_Isotope = QVariant(ls_Key.split("-").last()).toInt();
                        lk_MatchValues << tk_DoublePair(mk_UnlabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey][li_Isotope].first, lr_Peak.md_PeakIntensity);
                    }
                    double ld_Factor = leastSquaresFit(lk_MatchValues);
                    lr_ScanResult.md_UnlabeledProfileScale = ld_Factor;
                    
                    double ld_Error = 0.0;
                    foreach (tk_DoublePair lk_Pair, lk_MatchValues)
                    {
                        double ld_EnvelopeHeight = lk_Pair.first;
                        double ld_PeakHeight = lk_Pair.second;
                        ld_PeakHeight /= ld_Factor;
                        ld_Error += pow(ld_PeakHeight - ld_EnvelopeHeight, 2.0);
                    }
                    lr_ScanResult.md_UnlabeledError = ld_Error;
                    if (mb_UseIsotopeEnvelopes)
                        lr_ScanResult.md_AmountUnlabeled = ld_Factor;
                }
                
                if (!lk_LabeledKeys.empty())
                {
                    QList<tk_DoublePair> lk_MatchValues;
                    foreach (QString ls_Key, lk_LabeledKeys)
                    {
                        int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
                        r_Peak lr_Peak = lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];
//                         if (me_AmountEstimation == r_AmountEstimation::Area)
//                             lr_ScanResult.md_AmountLabeled += lr_Peak.md_PeakArea;
//                         else if (me_AmountEstimation == r_AmountEstimation::Intensity)
                        // :TODO: now this always uses peak height, make this an option (area/intensity)
                        lr_ScanResult.md_AmountLabeled += lr_Peak.md_PeakIntensity;
                        lr_ScanResult.md_Snr = std::min<double>(lr_Peak.md_Snr, lr_ScanResult.md_Snr);
                        lr_ScanResult.mk_LabeledPeaks.append(lr_Peak);
                        int li_Isotope = QVariant(ls_Key.split("-").last()).toInt();
                        lk_MatchValues << tk_DoublePair(mk_LabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey][li_Isotope].first, lr_Peak.md_PeakIntensity);
                    }
                    double ld_Factor = leastSquaresFit(lk_MatchValues);
                    lr_ScanResult.md_LabeledProfileScale = ld_Factor;
                    
                    double ld_Error = 0.0;
                    foreach (tk_DoublePair lk_Pair, lk_MatchValues)
                    {
                        double ld_EnvelopeHeight = lk_Pair.first;
                        double ld_PeakHeight = lk_Pair.second;
                        ld_PeakHeight /= ld_Factor;
                        ld_Error += pow(ld_PeakHeight - ld_EnvelopeHeight, 2.0);
                    }
                    lr_ScanResult.md_LabeledError = ld_Error;
                    if (mb_UseIsotopeEnvelopes)
                        lr_ScanResult.md_AmountLabeled = ld_Factor;
                }
                
                // check whether QE is good
                bool lb_GoodQE = true;
                if (lr_ScanResult.md_Snr < md_MinSnr)
                    lb_GoodQE = false;
                
                if (lb_GoodQE)
                {
                    if (mk_CsvOutStream.device())
                    {
                        mk_CsvOutStream << "\"" << ms_CurrentSpot << "\""
                            << ",\"" << ar_Scan.ms_Id << "\""
                            << ",\"" << ls_Peptide << "\""
                            << "," << lr_ScanResult.md_AmountUnlabeled
                            << "," << lr_ScanResult.md_AmountLabeled
                            << "," << ar_Scan.md_RetentionTime
                            << "," << lr_ScanResult.mi_Charge
                            << ",\"" << ar_Scan.ms_FilterLine << "\""
                            << "," << lr_ScanResult.md_Snr
                            << endl;
                    }
                
                    if (mk_XhtmlOutStream.device())
                    {
                        mk_XhtmlOutStream << QString("\n<!-- BEGIN PEPTIDE %1 -->\n").arg(ls_Peptide); 
                        mk_XhtmlOutStream << "<tr>"
                            << "<td>" << ms_CurrentSpot << "</td>"
                            << "<td>" << ar_Scan.ms_Id << "</td>"
                            << "<td>" << ls_Peptide << "</td>"
                            << "<td>" << lr_ScanResult.md_AmountUnlabeled << "</td>"
                            << "<td>" << lr_ScanResult.md_AmountLabeled << "</td>"
                            << "<td>" << ar_Scan.md_RetentionTime << "</td>"
                            << "<td>" << lr_ScanResult.mi_Charge << "</td>"
                            << "<td>" << ar_Scan.ms_FilterLine << "</td>"
                            << "<td>" << lr_ScanResult.md_Snr << "</td>"
                            << "</tr>"
                            << endl;
                        QString ls_Svg = this->renderScanAsSvg(ar_Scan, lr_ScanResult);
                        ls_Svg.remove(QRegExp("<\\?xml.+\\?>"));
                        ls_Svg.replace(QRegExp("width=\\\"[^\\\"]*\\\"\\s+height=\\\"[^\\\"]*\\\""), "width='950' height='238'");
                        mk_XhtmlOutStream << "<div style='background-color: #fff;' width='950' height='238'>";
                        mk_XhtmlOutStream << ls_Svg;
                        mk_XhtmlOutStream << "</div>" << endl;
                        mk_XhtmlOutStream << QString("\n<!-- END PEPTIDE %1 -->\n").arg(ls_Peptide); 
                    }
                }
            }
        }
	}
    */
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
