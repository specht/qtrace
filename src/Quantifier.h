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

#pragma once
#include <QtCore>
#include <ptb/ScanIterator.h>


#define HYDROGEN 1.007947
#define WATER 18.01057
//#define NEUTRON 1.008664915
// TODO: verify this neutron mass
#define NEUTRON 1.002
#define HEAVY_ARGININE 6.02013
#define HEAVY_PROLINE 5.016775
#define HEAVY_NITROGEN 1.003355

struct r_LabelType
{
	enum Enumeration {
		HeavyArginine = 0,
		HeavyArginineAndProline = 1,
		N15Labeling = 2,
		Size
	};
};

struct r_QuantitationFailureReason
{
	enum Enumeration {
		NoMatchedTargetMass,
		IsotopePeaksMissing,
		ForbiddenPeakPresent,
		LowSnr,
		Success,
		Size
	};
};


struct r_ScanQuantitationResult
{
	bool mb_IsGood;
	double md_AmountUnlabeled;
	double md_AmountLabeled;
	double md_Ratio;
	double md_Snr;
	double md_RetentionTime;
	QList<double> mk_UnlabeledTargetMz;
	QList<double> mk_LabeledTargetMz;
	double md_MinMz;
	double md_MaxMz;
	int mi_Charge;
	QString ms_ScanId;
	QList<r_Peak> mk_UnlabeledPeaks;
	QList<r_Peak> mk_LabeledPeaks;
	
	// scan data
	QString ms_ScanHashKey;
};


struct r_Bucket
{
	r_Bucket(int ai_Start, int ai_Length)
		: mi_Start(ai_Start)
		, mi_Length(ai_Length)
	{
	}
	
	int mi_Start;
	int mi_Length;
	QList<int> mk_Entries;
};
		

// peptide => scan results
typedef QHash<QString, QList<r_ScanQuantitationResult> > tk_SpotResults;


class k_Quantifier: public k_ScanIterator
{
public:
	k_Quantifier(r_LabelType::Enumeration ae_LabelType = r_LabelType::HeavyArginineAndProline,
				 r_ScanType::Enumeration ae_ScanType = r_ScanType::All,
				 QList<tk_IntPair> ak_MsLevels = QList<tk_IntPair>() << tk_IntPair(0, 0x10000),
				 int ai_IsotopeCount = 3, int ai_MinCharge = 2, int ai_MaxCharge = 3, 
				 double ad_MinSnr = 2.0, double ad_MassAccuracy = 5.0,
				 double ad_ExcludeMassAccuracy = 10.0,
				 QIODevice* ak_CsvOutDevice_ = NULL, QIODevice* ak_XhtmlOutDevice_ = NULL,
				 bool ab_PrintStatistics = false);
	virtual ~k_Quantifier();
	
	// quantify takes a list of spectra files and a hash of (peptide => protein) entries
	virtual void quantify(QStringList ak_SpectraFiles, QStringList ak_Peptides);
	virtual void handleScan(r_Scan& ar_Scan);
	virtual void progressFunction(QString as_ScanId, bool ab_InterestingScan);
	
	virtual QString renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult);
	
protected:
	virtual double calculatePeptideMass(QString as_Peptide, int ai_Charge);
	virtual double scale(const double ad_Value) const;
	virtual void calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_);
	virtual r_ScanQuantitationResult 
		checkResult(QHash<int, r_Peak> ak_LightPeaksInclude, 
					 QHash<int, r_Peak> ak_HeavyPeaksInclude,
					 QHash<int, r_Peak> ak_LightPeaksExclude, 
					 QHash<int, r_Peak> ak_HeavyPeaksExclude,
					 r_Scan& ar_Scan, QString as_Peptide, 
					 int ai_Charge, 
					 QList<double> ak_UnlabeledTargetMz, 
					 QList<double> ak_LabeledTargetMz);
	void fitGaussian(double* a_, double* b_, double* c_, double x0, double y0, 
					 double x1, double y1, double x2, double y2);
	double gaussian(double x, double a, double b, double c);
	void updateFailureReason(const QString& as_Key, r_QuantitationFailureReason::Enumeration ae_Value);

/*
	peptide:
		spectra file/charge:
			ratio:
			certainty:
			more info:
*/
	int mi_WatchIsotopesCount;
	QTextStream mk_CsvOutStream;
	QTextStream mk_XhtmlOutStream;
	r_LabelType::Enumeration me_LabelType;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_MassAccuracy;
	double md_ExcludeMassAccuracy;
	double md_ElutionProfilePeakWidth;
	bool mb_PrintStatistics;
	QList<double> mk_AllTargetMasses;
	QString ms_CurrentSpot;
	QStringList mk_Peptides;
	QHash<char, double> mk_AminoAcidWeight;
	QHash<char, int> mk_AminoAcidNitrogenCount;
	unsigned int mui_QuantitationResultCount;
	
	// peptide-charge-label-isotope
	QHash<QString, int> mk_TargetMzIndex;
	
	QHash<QString, r_QuantitationFailureReason::Enumeration> mk_ScanFailureReason;
	QHash<QString, int> mk_LabeledEnvelopeCountForPeptide;
	
	//QHash<QString, QList<QList<double> > > mk_ElutionProfile;
};
