/*
Copyright (c) 2007-2009 Michael Specht

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


k_Quantifier::k_Quantifier(QString as_Label,
						   r_ScanType::Enumeration ae_ScanType,
                           r_AmountEstimation::Enumeration ae_AmountEstimation,
						   QList<tk_IntPair> ak_MsLevels,
						   int ai_MinCharge, int ai_MaxCharge, 
						   double ad_MinSnr, double ad_MassAccuracy,
                           double ad_RequireAbundance,
                           double ad_ConsiderAbundance,
                           double ad_MaxFitError,
						   QIODevice* ak_CsvOutDevice_, QIODevice* ak_XhtmlOutDevice_,
						   bool ab_CheckForbiddenPeak,
						   bool ab_PrintStatusMessages, bool ab_LogScale)
	: k_ScanIterator(ae_ScanType, ak_MsLevels)
	, ms_Label(as_Label)
    , me_AmountEstimation(ae_AmountEstimation)
	, mk_CsvOutStream(ak_CsvOutDevice_)
	, mk_XhtmlOutStream(ak_XhtmlOutDevice_)
	, mi_MinCharge(ai_MinCharge)
	, mi_MaxCharge(ai_MaxCharge)
	, md_MinSnr(ad_MinSnr)
	, md_MassAccuracy(ad_MassAccuracy)
    , md_RequireAbundance(ad_RequireAbundance)
    , md_ConsiderAbundance(ad_ConsiderAbundance)
    , md_MaxFitError(ad_MaxFitError)
	, mb_CheckForbiddenPeak(ab_CheckForbiddenPeak)
	, mb_PrintStatusMessages(ab_PrintStatusMessages)
    , mb_LogScale(ab_LogScale)
{
    md_HydrogenMass = mk_IsotopeEnvelope.mk_BaseIsotopeMass["H"];
    double ld_OxygenMass = mk_IsotopeEnvelope.mk_BaseIsotopeMass["O"];
    md_WaterMass = md_HydrogenMass * 2.0 + ld_OxygenMass;
    
	// test parameter sanity
	if (mi_MinCharge < 1)
	{
		printf("Error: minimum charge (%d) out of range.\n", mi_MinCharge);
		exit(1);
	}
	if (mi_MaxCharge < 1)
	{
		printf("Error: maximum charge (%d) out of range.\n", mi_MaxCharge);
		exit(1);
	}
	if (mi_MaxCharge < mi_MinCharge)
	{
		printf("Error: minimum charge (%d) bigger than maximum charge(%d).\n", mi_MinCharge, mi_MaxCharge);
		exit(1);
	}
	if (md_MassAccuracy < 0.0)
	{
		printf("Error: mass accuracy must not be less than zero.\n");
		exit(1);
	}
    if (md_RequireAbundance < 0.0 || md_RequireAbundance > 1.0)
    {
        printf("Error: require abundance out of range (0.0 - 1.0).\n");
        exit(1);
    }
    if (md_ConsiderAbundance < 0.0 || md_ConsiderAbundance > 1.0)
    {
        printf("Error: consider abundance out of range (0.0 - 1.0).\n");
        exit(1);
    }
    if (md_ConsiderAbundance > md_RequireAbundance)
    {
        printf("Error: consider abundance must not be greater than require abundance.\n");
        exit(1);
    }
	
	Q_INIT_RESOURCE(qtrace);
    
	QFile lk_File(":ext/proteomics-knowledge-base/amino-acids.csv");
	if (!lk_File.open(QFile::ReadOnly))
    {
        fprintf(stderr, "Error: Unable to open amino-acids.csv.\n");
        exit(1);
    }
	QTextStream lk_TextStream(&lk_File);

    QString ls_Header = lk_TextStream.readLine().toLower().trimmed();
    QStringList lk_Header = ls_Header.split(",");
    int li_LetterIndex = lk_Header.indexOf("single letter code");
    int li_CompositionIndex = lk_Header.indexOf("composition");
    int li_MassIndex = lk_Header.indexOf("monoisotopic mass");
    if (li_LetterIndex < 0 || li_CompositionIndex < 0 || li_MassIndex < 0)
    {
        fprintf(stderr, "Error: Something is wrong with amino-acids.csv.\n");
        exit(1);
    }
	while (!lk_TextStream.atEnd())
	{
		QString ls_Line = lk_TextStream.readLine().trimmed();
		if (ls_Line.isEmpty())
			continue;

		QStringList lk_List;
		foreach (QString ls_Entry, ls_Line.split(QChar(',')))
		{
			if (ls_Entry.at(0) == QChar('"') && ls_Entry.at(ls_Entry.length() - 1) == QChar('"'))
				ls_Entry = ls_Entry.mid(1, ls_Entry.length() - 2);
			lk_List << ls_Entry;
		}
		char lc_AminoAcid = lk_List[li_LetterIndex][0].toAscii();
		mk_AminoAcidWeight[lc_AminoAcid] = lk_List[li_MassIndex].toDouble();
		QString ls_Composition = lk_List[li_CompositionIndex];
		ls_Composition.remove("\"");
		QHash<QString, int> lk_ElementCount;
		while (!ls_Composition.isEmpty())
		{
			char lc_Element = ls_Composition.at(0).toAscii();
			lk_ElementCount[QString(lc_Element)] = 1;
			ls_Composition = ls_Composition.right(ls_Composition.length() - 1);
			if (ls_Composition.isEmpty())
				break;
			QString ls_Number;
			char lc_NextChar = ls_Composition.at(0).toAscii();
			while (lc_NextChar >= '0' && lc_NextChar <= '9')
			{
				ls_Number += QChar(lc_NextChar);
				ls_Composition = ls_Composition.right(ls_Composition.length() - 1);
				if (ls_Composition.isEmpty())
					break;
				lc_NextChar = ls_Composition.at(0).toAscii();
			}
			if (!ls_Number.isEmpty())
				lk_ElementCount[QString(lc_Element)] = ls_Number.toInt();
		}
        mk_AminoAcidComposition[lc_AminoAcid] = lk_ElementCount;
	}
	lk_File.close();
    
    parseLabel();
}


k_Quantifier::~k_Quantifier()
{
}


bool compareScanQuantitationBySnr(const r_ScanQuantitationResult& a, const r_ScanQuantitationResult& b)
{
	return a.md_Snr > b.md_Snr;
}


bool compareScanQuantitationByRetentionTime(const r_ScanQuantitationResult& a, const r_ScanQuantitationResult& b)
{
	return a.md_RetentionTime < b.md_RetentionTime;
}


bool sortByMz(const QPair<double, QString>& a, const QPair<double, QString>& b)
{
	return a.first < b.first;
}


QString dtos(double ad_Value)
{
	static char lc_String_[2048];
	sprintf(lc_String_, "%1.2f", ad_Value);
	return QString(lc_String_);
}


void k_Quantifier::quantify(QStringList ak_SpectraFiles, QStringList ak_Peptides)
{
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

/*                        fprintf(stderr, "A+%d %c %9.5f / %9.5f\n", li_Isotope, 
                                lk_NormalizedEnvelope[li_Isotope].first >= REQUIRE_ABUNDANCE ? '*' : lk_NormalizedEnvelope[li_Isotope].first >= CONSIDER_ABUNDANCE ? '+' : ' ', 
                                lk_NormalizedEnvelope[li_Isotope].first, ld_Mz);*/
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

                        /*
                        fprintf(stderr, "A*+%d %c %9.5f / %9.5f\n", li_Isotope, 
                                lk_NormalizedEnvelopeN15[li_Isotope].first > REQUIRE_ABUNDANCE ? '*' : lk_NormalizedEnvelopeN15[li_Isotope].first > CONSIDER_ABUNDANCE ? '+' : ' ', 
                                lk_NormalizedEnvelopeN15[li_Isotope].first, ld_Mz);
                        */
                    }
                }
            }
        }
    }
    
	if (mb_PrintStatusMessages)
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
	
	qSort(lk_TempList.begin(), lk_TempList.end(), sortByMz);

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
		lk_File.close();
	}
		
	// parse all bands
	foreach (QString ls_Path, ak_SpectraFiles)
	{
		ms_CurrentSpot = QFileInfo(ls_Path).baseName();
		
		// parse spot
		if (!lk_PeptidesActuallySearchedFor.empty())
			this->parseFile(ls_Path);
		
		if (mb_PrintStatusMessages)
			printf(" done.\n");
	}

	if (mk_XhtmlOutStream.device())
	{
		QFile lk_File(":res/qtrace-xhtml-footer.xhtml.part");
		lk_File.open(QIODevice::ReadOnly);
		QByteArray lk_Content = lk_File.readAll();
		mk_XhtmlOutStream << QString(lk_Content);
		lk_File.close();
	}
}


void k_Quantifier::handleScan(r_Scan& ar_Scan)
{
/*	if (QVariant(ar_Scan.ms_Id).toInt() != 5117)
		return;*/
	
	if (ar_Scan.mr_Spectrum.mi_PeaksCount == 0)
	{
		printf("Warning: Empty spectrum (scan #%s @ %1.2f minutes)!\n", ar_Scan.ms_Id.toStdString().c_str(), ar_Scan.md_RetentionTime);
		return;
	}
	
	// find all peaks in this spectrum
	QList<r_Peak> lk_AllPeaks = k_ScanIterator::findAllPeaks(ar_Scan.mr_Spectrum);
// 	printf("all peaks: %d\n", lk_AllPeaks.size());
	
	// match all target m/z values simultaneously to this spectrum's peaks
	// create root bucket
	r_Bucket lr_RootBucket(0, lk_AllPeaks.size());
	for (int i = 0; i < mk_AllTargetMasses.size(); ++i)
		lr_RootBucket.mk_Entries.push_back(i);
	
	QList<r_Bucket> lk_Buckets;
	lk_Buckets.push_back(lr_RootBucket);
	
	// after parallel searching, this list will contain for 
	// each target m/z value the peak which is closest to it
 	QHash<int, int> lk_PeakForTargetMz;
	
	// repeat until no more buckets are left to be processed
	while (!lk_Buckets.empty())
	{
		// put all target m/z values into the appropriate bucket
		QList<r_Bucket> lk_NewBuckets;
		foreach (r_Bucket lr_Bucket, lk_Buckets)
		{
			if (lr_Bucket.mi_Length <= 1)
			{
				for (int i = 0; i < lr_Bucket.mk_Entries.size(); ++i)
				{
					/*
					printf("%1.6f is close to %1.6f\n", 
						mk_AllTargetMasses[lr_Bucket.mk_Entries[i]],
						ar_Scan.mr_Spectrum.md_MzValues_[lr_Bucket.mi_Start]);
					*/
					lk_PeakForTargetMz[lr_Bucket.mk_Entries[i]] = lr_Bucket.mi_Start;
				}
			}
			else
			{
				// split this bucket and create left and right children
				int li_HalfSize = lr_Bucket.mi_Length / 2;
				r_Bucket lr_LeftChild(lr_Bucket.mi_Start, li_HalfSize);
				r_Bucket lr_RightChild(lr_Bucket.mi_Start + li_HalfSize, lr_Bucket.mi_Length - li_HalfSize);
				double ld_LeftBorder = lk_AllPeaks[lr_RightChild.mi_Start - 1].md_PeakMz;
				double ld_RightBorder = lk_AllPeaks[lr_RightChild.mi_Start].md_PeakMz;
				double ld_Razor = (ld_LeftBorder + ld_RightBorder) * 0.5;
				
				// now determine which target m/z entries go into the 
				// left child and which go into the right child
				
				// TODO: this could be sped up because everything is sorted
				for (int i = 0; i < lr_Bucket.mk_Entries.size(); ++i)
				{
					if (mk_AllTargetMasses[lr_Bucket.mk_Entries[i]] > ld_Razor)
						lr_RightChild.mk_Entries.push_back(lr_Bucket.mk_Entries[i]);
					else
						lr_LeftChild.mk_Entries.push_back(lr_Bucket.mk_Entries[i]);
				}
				
				if (!lr_LeftChild.mk_Entries.empty())
					lk_NewBuckets.push_back(lr_LeftChild);
				if (!lr_RightChild.mk_Entries.empty())
					lk_NewBuckets.push_back(lr_RightChild);
			}
		}
		lk_Buckets = lk_NewBuckets;
	}
	
/*	foreach (QString s, mk_TargetMzIndex.keys())
	{
		double ld_ScanMz = lk_AllPeaks[lk_PeakForTargetMz[mk_TargetMzIndex[s]]].md_PeakMz;
		double ld_TargetMz = mk_AllTargetMasses[mk_TargetMzIndex[s]];
		double ld_Error = fabs(ld_ScanMz - ld_TargetMz) / ld_TargetMz * 1000000.0;
		double ld_Snr = lk_AllPeaks[lk_PeakForTargetMz[mk_TargetMzIndex[s]]].md_Snr;
		printf("%s: %1.6f - %1.6f (%1.2f ppm, SNR %1.2f)\n", s.toStdString().c_str(), ld_ScanMz, ld_TargetMz, ld_Error, ld_Snr);
	}*/
	
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
/*                printf("%s %d+ %d %d\n", ls_Peptide.toStdString().c_str(),
                       li_Charge, lb_FoundLightEnvelope, lb_FoundHeavyEnvelope);*/
                
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
                        if (me_AmountEstimation == r_AmountEstimation::Area)
                            lr_ScanResult.md_AmountUnlabeled += lr_Peak.md_PeakArea;
                        else if (me_AmountEstimation == r_AmountEstimation::Intensity)
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
                    if (me_AmountEstimation == r_AmountEstimation::Profile)
                        lr_ScanResult.md_AmountUnlabeled = ld_Factor;
                }
                
                if (!lk_LabeledKeys.empty())
                {
                    QList<tk_DoublePair> lk_MatchValues;
                    foreach (QString ls_Key, lk_LabeledKeys)
                    {
                        int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
                        r_Peak lr_Peak = lk_AllPeaks[lk_PeakForTargetMz[li_TargetMzIndex]];
                        if (me_AmountEstimation == r_AmountEstimation::Area)
                            lr_ScanResult.md_AmountLabeled += lr_Peak.md_PeakArea;
                        else if (me_AmountEstimation == r_AmountEstimation::Intensity)
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
                    if (me_AmountEstimation == r_AmountEstimation::Profile)
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
    return;
	
	// now check for each peptide/charge whether all peaks are there
	foreach (QString ls_Peptide, mk_Peptides)
	{
		for (int li_Charge = mi_MinCharge; li_Charge <= mi_MaxCharge; ++li_Charge)
		{
			QHash<int, r_Peak> lk_LightPeaksInclude;
			QHash<int, r_Peak> lk_HeavyPeaksInclude;
			QHash<int, r_Peak> lk_LightPeaksExclude;
			QHash<int, r_Peak> lk_HeavyPeaksExclude;
			QList<double> lk_UnlabeledTargetMz;
			QList<double> lk_LabeledTargetMz;
			
			QList<double> lk_ForbiddenMz;
			
/*			for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
			{
				QString ls_Key = QString("%1-%2-unlabeled-%3").
					arg(ls_Peptide).arg(li_Charge).arg(li_Isotope);
				int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
				lk_UnlabeledTargetMz.push_back(mk_AllTargetMasses[li_TargetMzIndex]);
				if (lk_IncludePeakForTargetMz.contains(li_TargetMzIndex))
					lk_LightPeaksInclude[li_Isotope] = lk_AllPeaks[lk_IncludePeakForTargetMz[li_TargetMzIndex]];
				if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
					lk_LightPeaksExclude[li_Isotope] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
			} */
/*			for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[ls_Peptide]; ++k)
			{
				for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
				{
					QString ls_Key = QString("%1-%2-labeled-%3-%4").
						arg(ls_Peptide).arg(li_Charge).arg(k).arg(li_Isotope);
					int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
					lk_LabeledTargetMz.push_back(mk_AllTargetMasses[li_TargetMzIndex]);
					if (lk_IncludePeakForTargetMz.contains(li_TargetMzIndex))
						lk_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope] = lk_AllPeaks[lk_IncludePeakForTargetMz[li_TargetMzIndex]];
					if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
						lk_HeavyPeaksExclude[k * mi_WatchIsotopesCount + li_Isotope] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
				}
				QString ls_Key = QString("%1-%2-forbidden-labeled-%3").arg(ls_Peptide).arg(li_Charge).arg(k);
				int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
				if (mb_CheckHeavyForbiddenPeaks)
					lk_ForbiddenMz << mk_AllTargetMasses[li_TargetMzIndex];
				if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
					lk_HeavyPeaksExclude[-1 - k] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
			}*/
			QString ls_Key = QString("%1-%2-forbidden").arg(ls_Peptide).arg(li_Charge);
			int li_TargetMzIndex = mk_TargetMzIndex[ls_Key];
			if (mb_CheckForbiddenPeak)
				lk_ForbiddenMz << mk_AllTargetMasses[li_TargetMzIndex];
			if (lk_ExcludePeakForTargetMz.contains(li_TargetMzIndex))
				lk_LightPeaksExclude[-1] = lk_AllPeaks[lk_ExcludePeakForTargetMz[li_TargetMzIndex]];
			
			QList<double> lk_TargetMz = lk_UnlabeledTargetMz + lk_LabeledTargetMz;
			
			r_ScanQuantitationResult lr_ScanResult = 
				this->checkResult(lk_LightPeaksInclude, lk_HeavyPeaksInclude,
								   lk_LightPeaksExclude, lk_HeavyPeaksExclude,
								   ar_Scan, ls_Peptide, li_Charge, 
								   lk_TargetMz, lk_ForbiddenMz);
			if (lr_ScanResult.mb_IsGood)
			{
				lr_ScanResult.ms_ScanHashKey = QString("%1.%2.%3.%4").arg(ms_CurrentSpot).arg(ar_Scan.ms_Id).arg(ls_Peptide).arg(li_Charge);
				lr_ScanResult.md_RetentionTime = ar_Scan.md_RetentionTime;
                lr_ScanResult.ms_Peptide = ls_Peptide;
				lr_ScanResult.mi_Charge = li_Charge;
				lr_ScanResult.ms_ScanId = ar_Scan.ms_Id;
				/*
				mk_ScanHash.insert(lr_ScanResult.ms_ScanHashKey, r_Scan(ar_Scan));
				if (!mk_SpotResults.contains(ls_Peptide))
					mk_SpotResults[ls_Peptide] = QList<r_ScanQuantitationResult>();
				mk_SpotResults[ls_Peptide].push_back(lr_ScanResult);

				*/
			
				if (mk_CsvOutStream.device())
				{
// 						mk_CsvOutStream << "id,filename,scan id,peptide,amount light,amount heavy,retention time,charge,filter line,snr" << endl;
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


void k_Quantifier::progressFunction(QString as_ScanId, bool)
{
	if (mb_PrintStatusMessages)
		printf("\r%s: scan #%s...", ms_CurrentSpot.toStdString().c_str(), as_ScanId.toStdString().c_str());
}


double k_Quantifier::calculatePeptideMass(QString as_Peptide, int ai_Charge)
{
	double ld_Mass = md_WaterMass;
	for (int i = 0; i < as_Peptide.length(); ++i)
		ld_Mass += mk_AminoAcidWeight[as_Peptide.at(i).toAscii()];
	// UNSURE ABOUT THIS VALUE BUT OK WITH BIANCAS EXAMPLES, should be 1.0078250
	return (ld_Mass + md_HydrogenMass * ai_Charge) / ai_Charge;
}


QString k_Quantifier::renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult)
{
	double ld_Ratio = 4.0;
	double ld_Width = 950.0;
	double ld_Height = ld_Width / ld_Ratio;
	double ld_BorderTop = 4.0;
	double ld_BorderRight = 16.0;
	double ld_BorderBottom = 4.0;
	double ld_BorderLeft = 16.0;
	
	if (ar_QuantitationResult.md_MinMz == 0.0)
		ar_QuantitationResult.md_MinMz = ar_Scan.mr_Spectrum.md_MzValues_[0];
	if (ar_QuantitationResult.md_MaxMz == 0.0)
		ar_QuantitationResult.md_MaxMz = ar_Scan.mr_Spectrum.md_MzValues_[ar_Scan.mr_Spectrum.mi_PeaksCount - 1];
		
	int li_Start = 0;
	int li_End = ar_Scan.mr_Spectrum.mi_PeaksCount - 1;
	
	// adjust start and end values
	while (li_Start < ar_Scan.mr_Spectrum.mi_PeaksCount - 1 && ar_Scan.mr_Spectrum.md_MzValues_[li_Start] < ar_QuantitationResult.md_MinMz)
		++li_Start;
	while (li_End > 0 && ar_Scan.mr_Spectrum.md_MzValues_[li_End] > ar_QuantitationResult.md_MaxMz)
		--li_End;
	
	if (li_Start > 0)
		--li_Start;
	if (li_End < ar_Scan.mr_Spectrum.mi_PeaksCount - 1)
		++li_End;
	if (li_End < ar_Scan.mr_Spectrum.mi_PeaksCount - 1)
		++li_End;
		
	if (li_Start > li_End)
	{
		// oops! start and end got crossed.
		li_Start = 0;
		li_End = 0;
	}
	
	// determine x and y extensions
	double xmin = ar_Scan.mr_Spectrum.md_MzValues_[li_Start];
	double xmax = ar_Scan.mr_Spectrum.md_MzValues_[li_End];
	xmin = ar_QuantitationResult.md_MinMz;
	xmax = ar_QuantitationResult.md_MaxMz;
	double ymax = ar_Scan.mr_Spectrum.md_IntensityValues_[li_Start];
	for (int i = li_Start + 1; i <= li_End; ++i)
		ymax = std::max<double>(ymax, ar_Scan.mr_Spectrum.md_IntensityValues_[i]);
	// ymax remains unscaled!
	double ymaxScaled = this->scale(1.0);
    
    double ymaxUnlabeled = 0.0;
    foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_UnlabeledPeaks)
        ymaxUnlabeled = std::max<double>(ymaxUnlabeled, lr_Peak.md_PeakIntensity / ymax);
	
	QBuffer lk_Buffer;
	lk_Buffer.open(QBuffer::WriteOnly);

    QList<double> lk_Lines;
    if (mb_LogScale)
        lk_Lines << 0.01 << 0.1 << 1.0;
    else
        lk_Lines << 0.2 << 0.4 << 0.6 << 0.8 << 1.0;
            
	double x0 = ld_BorderLeft;
	double y0 = ld_Height - ld_BorderBottom;
	double dx = (ld_Width - ld_BorderLeft - ld_BorderRight) / (xmax - xmin);
	double dy = -(ld_Height - ld_BorderTop - ld_BorderBottom);
    QString ls_Labels;
	
	{
		QSvgGenerator lk_Generator;
		lk_Generator.setSize(QSize((int)ld_Width, (int)ld_Height));
		lk_Generator.setOutputDevice(&lk_Buffer);
		
		QPainter lk_Painter(&lk_Generator);
		
		// fill background
		lk_Painter.fillRect(QRectF(0.0, 0.0, ld_Width, ld_Height), QBrush(Qt::white));
		
		// draw frame
		QPen lk_Pen;
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::BevelJoin);
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Painter.setPen(lk_Pen);
		
		lk_Painter.drawLine(QPointF(0.0, y0), QPointF(ld_Width, y0));
		//lk_Painter.drawLine(QPointF(0.0, y0), QPointF(ld_Width, y0));
		//lk_Painter.drawLine(QPointF(x0, y0), QPointF(x0, y0 + (ymax - ymin) * dy));
        
        QString ls_Peptide = ar_QuantitationResult.ms_Peptide;
        int li_Charge = ar_QuantitationResult.mi_Charge;
        QString ls_PeptideChargeKey = QString("%1-%2").arg(ls_Peptide).arg(li_Charge);

        tk_IsotopeEnvelope lk_IsotopeEnvelope;
        bool lb_DrawnOnePoint;
        double ld_LastMz;
        
        lk_IsotopeEnvelope = mk_UnlabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey];
        double ld_PeptideChargeBaseMz = calculatePeptideMass(ls_Peptide, li_Charge);
        lb_DrawnOnePoint = false;
        ld_LastMz = 0.0;
        QPainterPath lk_Path;
        for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
        {
            double ld_Abundance = lk_IsotopeEnvelope[i].first * ar_QuantitationResult.md_UnlabeledProfileScale;
            double ld_Mz = lk_IsotopeEnvelope[i].second / li_Charge;
            ld_Mz += ld_PeptideChargeBaseMz;
            double x = (ld_Mz - xmin) * dx + x0;
            double y = this->scale(ld_Abundance / ymax) / ymaxScaled;
            if (y > 0.001)
            {
                y = y * dy + y0;
                if (!lb_DrawnOnePoint)
                {
                    lk_Path.moveTo(QPointF(((ld_Mz - 0.5 / li_Charge) - xmin) * dx + x0, y0));
                    lb_DrawnOnePoint = true;
                    //lk_Path.quadTo(QPointF(((ld_Mz / li_Charge) - xmin) * dx + x0, y0), QPointF(x, y));
                }
                else
                {
                    //lk_Path.lineTo(QPointF(x, y));
                }
                lk_Path.lineTo(QPointF(x, y));
                ld_LastMz = ld_Mz;
            }
        }
        lk_Path.lineTo(QPointF(((ld_LastMz + 0.5 / li_Charge) - xmin) * dx + x0, y0));
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::BevelJoin);
        lk_Pen.setColor(QColor(128, 128, 128));
        lk_Painter.setPen(lk_Pen);
        lk_Painter.setBrush(QBrush(QColor(224, 224, 224, 128)));
        
        lk_Painter.drawPath(lk_Path);
		
        lk_IsotopeEnvelope = mk_LabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey];
        lb_DrawnOnePoint = false;
        ld_LastMz = 0.0;
        lk_Path = QPainterPath();
        for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
        {
            double ld_Abundance = lk_IsotopeEnvelope[i].first * ar_QuantitationResult.md_LabeledProfileScale;
            double ld_Mz = lk_IsotopeEnvelope[i].second / li_Charge;
            ld_Mz += ld_PeptideChargeBaseMz;
            double x = (ld_Mz - xmin) * dx + x0;
            double y = this->scale(ld_Abundance / ymax) / ymaxScaled;
            if (y > 0.001)
            {
                y = y * dy + y0;
                if (!lb_DrawnOnePoint)
                {
                    lk_Path.moveTo(QPointF(((ld_Mz - 0.5 / li_Charge) - xmin) * dx + x0, y0));
                    lb_DrawnOnePoint = true;
                    //lk_Path.quadTo(QPointF(((ld_Mz / li_Charge) - xmin) * dx + x0, y0), QPointF(x, y));
                }
                else
                {
                    //lk_Path.lineTo(QPointF(x, y));
                }
                lk_Path.lineTo(QPointF(x, y));
                ld_LastMz = ld_Mz;
            }
        }
        lk_Path.lineTo(QPointF(((ld_LastMz + 0.5 / li_Charge) - xmin) * dx + x0, y0));
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::BevelJoin);
        lk_Pen.setColor(QColor(128, 128, 128));
        lk_Painter.setPen(lk_Pen);
        lk_Painter.setBrush(QBrush(QColor(224, 224, 224, 128)));
        
        lk_Painter.drawPath(lk_Path);
        
        lk_Painter.setBrush(Qt::NoBrush);
		// draw spectrum
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::RoundJoin);
		lk_Pen.setColor(QColor(TANGO_CHAMELEON_2));
		lk_Pen.setStyle(Qt::DashLine);
		lk_Painter.setPen(lk_Pen);
		
		// draw target m/z lines
		foreach (double ld_Mz, ar_QuantitationResult.mk_TargetMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		
		lk_Pen.setColor(QColor(TANGO_SCARLET_RED_1));
		lk_Painter.setPen(lk_Pen);
		// draw forbidden m/z lines
		foreach (double ld_Mz, ar_QuantitationResult.mk_ForbiddenMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Painter.setPen(lk_Pen);
		
        foreach (double ld_Line, lk_Lines)
        {
            double y;
            y = this->scale(ld_Line) / ymaxScaled;
            y = y * dy + y0;
            lk_Painter.drawLine(QPointF(0.0, y), QPointF(ld_Width, y));
            ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana;' transform='rotate(90 %1 %2)'>%3%</text>\n").arg(0.0).arg(y + 3.0).arg((int)round(ld_Line * 100.0));
        }
		
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::RoundJoin);
		lk_Pen.setColor(QColor(128, 128, 128));
		lk_Pen.setStyle(Qt::SolidLine);
		lk_Painter.setPen(lk_Pen);
		
		QVector<QPointF> lk_Points;
		for (int i = li_Start; i <= li_End; ++i)
			lk_Points.append(QPointF((ar_Scan.mr_Spectrum.md_MzValues_[i] - xmin) * dx + x0, 
									this->scale(ar_Scan.mr_Spectrum.md_IntensityValues_[i] / ymax) / ymaxScaled * dy + y0));
		lk_Painter.drawPolyline(QPolygonF(lk_Points));
		
		foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_UnlabeledPeaks)
		{
			// draw Gaussian
			lk_Pen.setWidthF(1.0);
			lk_Pen.setJoinStyle(Qt::RoundJoin);
			lk_Pen.setColor(QColor(ar_QuantitationResult.md_UnlabeledError <= md_MaxFitError ? TANGO_SKY_BLUE_1 : TANGO_SCARLET_RED_1));
			lk_Pen.setStyle(Qt::SolidLine);
			lk_Painter.setPen(lk_Pen);
			
			QVector<QPointF> lk_Points;
			double w = sqrt(lr_Peak.md_GaussC) / 2.0;
			for (double x = lr_Peak.md_PeakMz - w; x <= lr_Peak.md_PeakMz + w; x += w * 0.05)
				lk_Points.append(QPointF((x - xmin) * dx + x0, 
										this->scale(gaussian(x, lr_Peak.md_GaussA, lr_Peak.md_GaussB, lr_Peak.md_GaussC) / ymax) / ymaxScaled * dy + y0));
			lk_Painter.drawPolyline(QPolygonF(lk_Points));
		}
        
        foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_LabeledPeaks)
        {
            // draw Gaussian
            lk_Pen.setWidthF(1.0);
            lk_Pen.setJoinStyle(Qt::RoundJoin);
            lk_Pen.setColor(QColor(ar_QuantitationResult.md_LabeledError <= md_MaxFitError ? TANGO_SKY_BLUE_1 : TANGO_SCARLET_RED_1));
            lk_Pen.setStyle(Qt::SolidLine);
            lk_Painter.setPen(lk_Pen);
            
            QVector<QPointF> lk_Points;
            double w = sqrt(lr_Peak.md_GaussC) / 2.0;
            for (double x = lr_Peak.md_PeakMz - w; x <= lr_Peak.md_PeakMz + w; x += w * 0.05)
                lk_Points.append(QPointF((x - xmin) * dx + x0, 
                                        this->scale(gaussian(x, lr_Peak.md_GaussA, lr_Peak.md_GaussB, lr_Peak.md_GaussC) / ymax) / ymaxScaled * dy + y0));
            lk_Painter.drawPolyline(QPolygonF(lk_Points));
        }
	}
	lk_Buffer.close();
	lk_Buffer.open(QBuffer::ReadOnly);
	lk_Buffer.seek(0);
	QString ls_Result = QString(lk_Buffer.readAll());
	lk_Buffer.close();
	foreach (r_Peak lr_Peak, (ar_QuantitationResult.mk_UnlabeledPeaks + ar_QuantitationResult.mk_LabeledPeaks))
	{
		double ld_X, ld_Y;
		ld_X = (lr_Peak.md_PeakMz - xmin) * dx + x0;
		ld_Y = this->scale(lr_Peak.md_PeakIntensity / ymax) / ymaxScaled * dy + y0;
		char lc_Number_[1024];
        lc_Number_[0] = 0;
        if (me_AmountEstimation == r_AmountEstimation::Area)
            sprintf(lc_Number_, "%.2g", lr_Peak.md_PeakArea);
        else if (me_AmountEstimation == r_AmountEstimation::Intensity)
            sprintf(lc_Number_, "%.2g", lr_Peak.md_PeakIntensity);
        if (strlen(lc_Number_) > 0)
            ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana; fill: #aaa;' transform='rotate(90 %1 %2)'>%3</text>\n").arg(ld_X + 3.0).arg(ld_Y + 3.0).arg(lc_Number_);
	}
	ls_Result.replace("</svg>", ls_Labels + "</svg>");
	return ls_Result;
}


double k_Quantifier::scale(const double ad_Value) const
{
    if (mb_LogScale)
        return log(ad_Value * 100.0 + 1.0);
    else
        return ad_Value;
}


void k_Quantifier::calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_)
{
	double ld_Mean = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_Mean += ld_Value;
		
	ld_Mean /= ak_Values.size();
	
	double ld_StandardDeviation = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_StandardDeviation += pow(ld_Value - ld_Mean, 2.0);
		
	ld_StandardDeviation /= ak_Values.size();
	ld_StandardDeviation = sqrt(ld_StandardDeviation);
	
	*ad_Mean_ = ld_Mean;
	*ad_StandardDeviation_ = ld_StandardDeviation;
}


r_ScanQuantitationResult 
k_Quantifier::checkResult(QHash<int, r_Peak> ak_LightPeaksInclude, 
						   QHash<int, r_Peak> ak_HeavyPeaksInclude,
						   QHash<int, r_Peak> ak_LightPeaksExclude, 
						   QHash<int, r_Peak> ak_HeavyPeaksExclude,
						   r_Scan& ar_Scan, QString as_Peptide, int ai_Charge,
						   QList<double> ak_TargetMz,
						   QList<double> ak_ForbiddenMz)
{
	r_ScanQuantitationResult lr_Result;
	lr_Result.mb_IsGood = false;
	
	// check whether all inclusion peaks have a good SNR
	QList<r_Peak> lk_Peaks = ak_LightPeaksInclude.values() + ak_HeavyPeaksInclude.values();
	
	// return if no inclusion peaks have been found
	if (lk_Peaks.empty())
		return lr_Result;
	
	double ld_MinSnr = lk_Peaks.first().md_Snr;
	for (int i = 1; i < lk_Peaks.size(); ++i)
		ld_MinSnr = std::min<double>(ld_MinSnr, lk_Peaks[i].md_Snr);
	
	if (ld_MinSnr < md_MinSnr)
		return lr_Result;
	
	lr_Result.md_Snr = ld_MinSnr;
	
	double ld_LightSum = 0.0;
	double ld_HeavySum = 0.0;
	
	int li_LightPeakIncludeCount = 0;
	int li_HeavyPeakIncludeCount = 0;
	int li_LightPeakExcludeCount = 0;
	int li_HeavyPeakExcludeCount = 0;
	
/*	for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
	{
		if (ak_LightPeaksInclude.contains(li_Isotope))
		{
            if (mb_UseArea)
                ld_LightSum += ak_LightPeaksInclude[li_Isotope].md_PeakArea;
            else
                ld_LightSum += ak_LightPeaksInclude[li_Isotope].md_PeakIntensity;
			++li_LightPeakIncludeCount;
		}
		if (ak_LightPeaksExclude.contains(li_Isotope))
			++li_LightPeakExcludeCount;*/
/*		for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
		{
			if (ak_HeavyPeaksInclude.contains(k * mi_WatchIsotopesCount + li_Isotope))
			{
                if (mb_UseArea)
                    ld_HeavySum += ak_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope].md_PeakArea;
                else
                    ld_HeavySum += ak_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope].md_PeakIntensity;
				++li_HeavyPeakIncludeCount;
			}
			if (ak_HeavyPeaksExclude.contains(k * mi_WatchIsotopesCount + li_Isotope))
				++li_HeavyPeakExcludeCount;
		}*/
// 	}
	
/*	if (li_LightPeakIncludeCount == mi_WatchIsotopesCount && 
		li_HeavyPeakIncludeCount == mi_WatchIsotopesCount * mk_LabeledEnvelopeCountForPeptide[as_Peptide])
	{
		// both isotope envelopes are complete, we get a ratio!
		if (mb_CheckLightForbiddenPeaks && ak_LightPeaksExclude.contains(-1))
			// don't quantify if forbidden peak is present!
			return lr_Result;
		if (mb_CheckHeavyForbiddenPeaks)
		{
			// and don't quantify if the heavy forbidden peak is present
			for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
				if (ak_HeavyPeaksExclude.contains(-1 - k))
					return lr_Result;
		}
		lr_Result.md_AmountUnlabeled = ld_LightSum;
		lr_Result.md_AmountLabeled = ld_HeavySum;
	}
	else if (li_LightPeakIncludeCount == mi_WatchIsotopesCount && 
			 li_HeavyPeakExcludeCount == 0)
	{
		// light state only!
		if (ak_LightPeaksExclude.contains(-1))
			return lr_Result;
		lr_Result.md_AmountUnlabeled = ld_LightSum;
		lr_Result.md_AmountLabeled = 0.0;
	}
	else if (li_LightPeakExcludeCount == 0 && 
			 li_HeavyPeakIncludeCount == mi_WatchIsotopesCount * mk_LabeledEnvelopeCountForPeptide[as_Peptide])
	{
		// heavy state only
		for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
			if (ak_HeavyPeaksExclude.contains(-1 - k))
				return lr_Result;
		lr_Result.md_AmountUnlabeled = 0.0;
		lr_Result.md_AmountLabeled = ld_HeavySum;
	}
	else
	{
		// no complete isotope envelope
		return lr_Result;
	}*/
	
	lr_Result.mk_TargetMz = ak_TargetMz;
	lr_Result.mk_ForbiddenMz = ak_ForbiddenMz;
	
	QList<double> lk_AllMz = ak_TargetMz + ak_ForbiddenMz;
	qSort(lk_AllMz);
	
	lr_Result.md_MinMz = lk_AllMz.first();
	lr_Result.md_MaxMz = lk_AllMz.last();
	
	lr_Result.mk_UnlabeledPeaks = ak_LightPeaksInclude.values();
	lr_Result.mk_LabeledPeaks = ak_HeavyPeaksInclude.values();
	lr_Result.mb_IsGood = true;
	
	return lr_Result;
}


void k_Quantifier::fitGaussian(double* a_, double* b_, double* c_, double x0, double y0, 
							   double x1, double y1, double x2, double y2)
{
	double A = 2.0 * (x1 - x0);
	double B = x0 * x0 - x1 * x1;
	double C = 2.0 * (x2 - x0);
	double D = x0 * x0 - x2 * x2;
	double lny0 = log(y0);
	double lny1 = log(y1);
	double lny2 = log(y2);
	double F = lny1 - lny0;
	double G = lny2 - lny0;
	double b = (F * D - B * G) / (A * G - F * C);
	double x0b = x0 - b;
	double x1b = x1 - b;
	double d = ((x0b * x0b) - (x1b * x1b)) / (lny1 - lny0);
	double a = y0 / (exp(-((x0b * x0b) / d)));
	double c = sqrt(d * 0.5);
	*a_ = a;
	*b_ = b;
	*c_ = c;
}
							   

double k_Quantifier::gaussian(double x, double a, double b, double c)
{
	return a * exp(-(pow(x - b, 2.0) / (2 * c * c)));
}


QHash<QString, int> k_Quantifier::compositionForPeptide(const QString& as_Peptide)
{
    QHash<QString, int> lk_Composition;
    lk_Composition["H"] = 2;
    lk_Composition["O"] = 1;
    for (int i = 0; i < as_Peptide.size(); ++i)
    {
        char lc_AminoAcid = as_Peptide.at(i).toAscii();
        QHash<QString, int> lk_AminoAcidComposition = mk_AminoAcidComposition[lc_AminoAcid];
        foreach (QString ls_Element, lk_AminoAcidComposition.keys())
        {
            if (!lk_Composition.contains(ls_Element))
                lk_Composition[ls_Element] = 0;
            lk_Composition[ls_Element] += lk_AminoAcidComposition[ls_Element];
        }
    }
    return lk_Composition;
}


double k_Quantifier::leastSquaresFit(QList<tk_DoublePair> ak_Pairs)
{
    double f = 0.0;
    double e = 0.0;
    foreach (tk_DoublePair lk_Pair, ak_Pairs)
    {
        f += lk_Pair.first * lk_Pair.second;
        e += lk_Pair.first * lk_Pair.first;
    }
    f /= e;
    return f;
}


void k_Quantifier::parseLabel()
{
    QStringList lk_Tokens = tokenize(ms_Label.toUpper());

    int li_State = 0;
    QSet<QString> lk_AminoAcidScope;
    
    // this hash contains special environmental conditions for every amino acid
    // if an amino acid is not contained, it comes from natural conditions!
    QHash<QString, tk_ArtificialEnvironment> lk_EnvironmentForAminoAcid;

    while (!lk_Tokens.empty())
    {
        QVariant::Type le_TokenType;
        QString ls_Token = fetchNextToken(&lk_Tokens, &le_TokenType);
        if (li_State == 0)
        {
            if (le_TokenType == QVariant::String)
            {
                // we found a scope-defining amino acid prefix
                lk_AminoAcidScope << ls_Token;
            }
            else
                li_State = 1;
        }
        if (li_State == 1)
        {
            if (le_TokenType != QVariant::Int)
            {
                printf("Error: Isotope number expected, but got %s.\n", ls_Token.toStdString().c_str());
                exit(1);
            }
            // fetch element
            if (lk_Tokens.empty())
            {
                printf("Error: Element expected, but string is at end.\n");
                exit(1);
            }
            
            int li_Isotope = ls_Token.toInt();
            
            ls_Token = fetchNextToken(&lk_Tokens, &le_TokenType);
            if (le_TokenType != QVariant::String)
            {
                printf("Error: Element expected, but got %s.\n", ls_Token.toStdString().c_str());
                exit(1);
            }
            
            QString ls_Element = ls_Token;
            double ld_Efficiency = 1.0;
            
            if (peekNextToken(lk_Tokens) == QVariant::Double)
            {
                // parse labeling efficiency
                ls_Token = fetchNextToken(&lk_Tokens, &le_TokenType);
                ld_Efficiency = ls_Token.toFloat();
                if (!(ls_Element == "C" || ls_Element == "N"))
                {
                    printf("Error: A labeling efficiency can only be specified for C and N.\n");
                    exit(1);
                }
            }
            
            if (lk_AminoAcidScope.empty() || lk_AminoAcidScope.contains("X"))
            {
                lk_AminoAcidScope = QSet<QString>();
                QString ls_AllAminoAcids = "GASPVTCLINDQKEMHFRYW";
                for (int i = 0; i < ls_AllAminoAcids.length(); ++i)
                    lk_AminoAcidScope << ls_AllAminoAcids.mid(i, 1);
            }
            foreach (QString ls_AminoAcid, lk_AminoAcidScope)
            {
                if (!lk_EnvironmentForAminoAcid.contains(ls_AminoAcid))
                    lk_EnvironmentForAminoAcid[ls_AminoAcid] = tk_ArtificialEnvironment();
                lk_EnvironmentForAminoAcid[ls_AminoAcid][ls_Element] = r_IsotopeAbundance(li_Isotope - mk_IsotopeEnvelope.mk_BaseIsotope[ls_Element], ld_Efficiency);
            }
            
            if (peekNextToken(lk_Tokens) != QVariant::Int)
            {
                lk_AminoAcidScope.clear();
                li_State = 0;
            }
        }
    }
    
    QMultiMap<QString, QString> lk_AminoAcidForDescription;
    
    foreach (QString ls_AminoAcid, lk_EnvironmentForAminoAcid.keys())
    {
//         printf("AE for %s:\n", ls_AminoAcid.toStdString().c_str());
        tk_ArtificialEnvironment lk_Environment = lk_EnvironmentForAminoAcid[ls_AminoAcid];
        QHash<QString, QList<double> > lk_Abundances;
        foreach (QString ls_Element, lk_Environment.keys())
        {
            r_IsotopeAbundance lr_IsotopeAbundance = lk_Environment[ls_Element];
//             printf("%d%s at %1.2f\n", lr_IsotopeAbundance.mi_Isotope, ls_Element.toStdString().c_str(), lr_IsotopeAbundance.mf_Efficiency);
            QList<double> lk_AbundanceNumbers;
            lk_AbundanceNumbers << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
            lk_AbundanceNumbers[lr_IsotopeAbundance.mi_Isotope] = lr_IsotopeAbundance.mf_Efficiency;
            if (lr_IsotopeAbundance.mf_Efficiency < 1.0)
                lk_AbundanceNumbers[1 - lr_IsotopeAbundance.mi_Isotope] = 1.0 - lr_IsotopeAbundance.mf_Efficiency;
            lk_Abundances[ls_Element] = lk_AbundanceNumbers;
        }
        mk_HeavyIsotopeEnvelopeForAminoAcid[ls_AminoAcid] = k_IsotopeEnvelope(lk_Abundances);
        
        // make description
        QStringList lk_ElementList = lk_Abundances.keys();
        qSort(lk_ElementList.begin(), lk_ElementList.end());
        QStringList lk_Description;
        foreach (QString ls_Element, lk_ElementList)
        {
            for (int i = 0; i < lk_Abundances[ls_Element].size(); ++i)
            {
                double ld_Abundance = lk_Abundances[ls_Element][i];
                if (ld_Abundance > 0.0)
                    lk_Description << QString("%1% %2%3")
                        .arg(ld_Abundance * 100, 0, 'f', 1)
                        .arg(i + mk_IsotopeEnvelope.mk_BaseIsotope[ls_Element])
                        .arg(ls_Element);
            }
        }
        lk_AminoAcidForDescription.insert(lk_Description.join(", "), ls_AminoAcid);
    }
    if (mb_PrintStatusMessages)
    {
        printf("Label composition:\n");
        foreach (QString ls_Description, lk_AminoAcidForDescription.uniqueKeys())
        {
            QStringList lk_AminoAcids = lk_AminoAcidForDescription.values(ls_Description);
            if (lk_AminoAcids.size() == 20)
                lk_AminoAcids = QStringList() << "all amino acids";
            printf("%s: %s\n", lk_AminoAcids.join(", ").toStdString().c_str(), ls_Description.toStdString().c_str());
        }
    }
}


QStringList k_Quantifier::tokenize(QString as_String)
{
    QStringList lk_Result;
    bool lb_Append = true;
    for (int i = 0; i < as_String.length(); ++i)
    {
        QChar lk_Char = as_String.at(i);
        if (lk_Char.isSpace())
            lb_Append = false;
        if (lk_Char.isDigit() || lk_Char == QChar('.'))
        {
            if (lk_Result.empty())
            {
                lk_Result << QString(lk_Char);
                lb_Append = true;
            }
            else
            {
                if (lb_Append & (lk_Result.last().at(0).isDigit() || lk_Result.last().at(0) == QChar('.')))
                    lk_Result.last() += lk_Char;
                else
                {
                    lk_Result << QString(lk_Char);
                    lb_Append = true;
                }
            }
        }
        else if (lk_Char.isLetter())
        {
            lk_Result << QString(lk_Char);
            lb_Append = true;
        }
    }
    return lk_Result;
}


QString k_Quantifier::fetchNextToken(QStringList* ak_StringList_, QVariant::Type* ae_Type_)
{
    QVariant::Type le_TokenType = QVariant::String;
    QString ls_Token = ak_StringList_->takeFirst();
    bool lb_Ok = false;
    ls_Token.toFloat(&lb_Ok);
    if (lb_Ok)
    {
        if (ls_Token.contains("."))
            le_TokenType = QVariant::Double;
        else
            le_TokenType = QVariant::Int;
    }
    *ae_Type_ = le_TokenType;
    return ls_Token;
}


QVariant::Type k_Quantifier::peekNextToken(QStringList ak_StringList)
{
    if (ak_StringList.empty())
        return QVariant::Invalid;
    QVariant::Type le_TokenType = QVariant::String;
    QString ls_Token = ak_StringList.first();
    bool lb_Ok = false;
    ls_Token.toFloat(&lb_Ok);
    if (lb_Ok)
    {
        if (ls_Token.contains("."))
            le_TokenType = QVariant::Double;
        else
            le_TokenType = QVariant::Int;
    }
    return le_TokenType;
}


tk_IsotopeEnvelope k_Quantifier::heavyEnvelopeForPeptide(QString as_Peptide)
{
    QHash<QString, int> lk_Composition = compositionForPeptide(as_Peptide);
    
    QHash<QString, QHash<QString, int> > lk_CompositionForAminoAcid;
    for (int i = 0; i < as_Peptide.length(); ++i)
    {
        QString ls_AminoAcid = as_Peptide.mid(i, 1);
        QString ls_Id = QString();
        if (mk_HeavyIsotopeEnvelopeForAminoAcid.contains(ls_AminoAcid))
            ls_Id = ls_AminoAcid;
        QHash<QString, int> lk_ThisComposition = compositionForPeptide(ls_AminoAcid);
        foreach (QString ls_Element, lk_ThisComposition.keys())
        {
            if (!lk_CompositionForAminoAcid[ls_Id].contains(ls_Element))
                lk_CompositionForAminoAcid[ls_Id][ls_Element] = 0;
            lk_CompositionForAminoAcid[ls_Id][ls_Element] += lk_ThisComposition[ls_Element];
        }
    }
    
    tk_IsotopeEnvelope lr_HeavyMass;
    bool lb_First = true;
    foreach (QString ls_Id, lk_CompositionForAminoAcid.keys())
    {
        k_IsotopeEnvelope* lk_UseEnvelope_ = &mk_IsotopeEnvelope;
        if (!ls_Id.isEmpty())
            lk_UseEnvelope_ = &mk_HeavyIsotopeEnvelopeForAminoAcid[ls_Id];
        if (lb_First)
            lr_HeavyMass = lk_UseEnvelope_->isotopeEnvelopeForComposition(lk_CompositionForAminoAcid[ls_Id]);
        else
            lr_HeavyMass = mk_IsotopeEnvelope.add(lr_HeavyMass, lk_UseEnvelope_->isotopeEnvelopeForComposition(lk_CompositionForAminoAcid[ls_Id]));
        lb_First = false;
    }
    
    return lr_HeavyMass;
}
