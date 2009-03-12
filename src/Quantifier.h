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
#include "ScanIterator.h"


#define HYDROGEN 1.007947
#define WATER 18.01057
//#define NEUTRON 1.008664915
// TODO: verify this neutron mass
#define NEUTRON 1.002
#define HEAVY_ARGININE 6.02013

struct r_QuantitationFailureReason
{
	enum Enumeration {
		NoMatchedTargetMass,
		IsotopePeaksMissing,
		UnlabeledIsotopePeaksNotDescending,
		ForbiddenPeakPresent,
		LowSnr,
		Success,
		Size
	};
};


struct r_Peak
{
	int mi_PeakIndex;
	int mi_LeftBorderIndex;
	int mi_RightBorderIndex;
	// attention, these are gauss peak m/z and intensity values! yay!
	double md_PeakMz;
	double md_PeakIntensity;
	double md_PeakArea;
	double md_OutsideBorderMaxIntensity;
	double md_Snr;
	double md_GaussA;
	double md_GaussB;
	double md_GaussC;
};


typedef QList<QPair<int, int> > tk_Highlight;


struct r_ScanQuantitationResult
{
	bool mb_IsGood;
	double md_AmountUnlabeled;
	double md_AmountLabeled;
	double md_Ratio;
	double md_Snr;
	double md_RetentionTime;
	tk_Highlight mk_Highlight;
	QList<double> mk_TargetMz;
	double md_MinMz;
	double md_MaxMz;
	int mi_Charge;
	QString ms_ScanId;
	QList<r_Peak> mk_Peaks;
	
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
	k_Quantifier(r_ScanType::Enumeration ae_ScanType = r_ScanType::All,
				 QList<tk_IntPair> ak_MsLevels = QList<tk_IntPair>() << tk_IntPair(0, 0x10000),
				 int ai_IsotopeCount = 3, int ai_MinCharge = 2, int ai_MaxCharge = 3, 
				 double ad_MinSnr = 2.0, double ad_MassAccuracy = 5.0,
				 QString as_SvgOutPath = QString(), QIODevice* ak_TextOutDevice_ = NULL,
				 QIODevice* ak_YamlOutDevice_ = NULL, bool ab_PrintStatistics = false);
	virtual ~k_Quantifier();
	
	// quantify takes a list of spectra files and a hash of (peptide => protein) entries
	virtual void quantify(QStringList ak_SpectraFiles, QStringList ak_Peptides);
	virtual void handleScan(r_Scan& ar_Scan);
	virtual void progressFunction(QString as_ScanId, bool ab_InterestingScan);
	
	virtual QString renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult);
	
protected:
	virtual double calculatePeptideMass(QString as_Peptide, int ai_Charge);
	virtual r_ScanQuantitationResult searchPeptide(r_Spectrum& ar_Spectrum, double ad_PeptideMz, int ai_Charge, double ad_ModificationMass);
	virtual QList<r_Peak> findPeaks(r_Spectrum& ar_Spectrum, QList<double> ak_TargetMasses);
	virtual QList<r_Peak> findAllPeaks(r_Spectrum& ar_Spectrum);
	virtual r_Peak findPeak(r_Spectrum& ar_Spectrum, int ai_Index);
	virtual double scale(const double ad_Value) const;
	virtual void calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_);
	virtual r_ScanQuantitationResult checkResult(QList<r_Peak> ak_Peaks, r_Peak* ar_ForbiddenPeak_, r_Scan& ar_Scan, QString as_Peptide, int ai_Charge, QList<double> ak_TargetMz);
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
	QTextStream mk_TextOutStream;
	QTextStream mk_YamlOutStream;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_MassAccuracy;
	double md_ElutionProfilePeakWidth;
	QString ms_SvgOutPath;
	bool mb_PrintStatistics;
	tk_SpotResults mk_SpotResults;
	QList<double> mk_AllTargetMasses;
	QString ms_CurrentSpot;
	QStringList mk_Peptides;
	QHash<char, double> mk_AminoAcidWeight;
	QHash<QString, r_Scan> mk_ScanHash;
	
	// peptide-charge-label-isotope
	QHash<QString, int> mk_TargetMzIndex;
	
	QHash<QString, r_QuantitationFailureReason::Enumeration> mk_ScanFailureReason;
	
	//QHash<QString, QList<QList<double> > > mk_ElutionProfile;
};
