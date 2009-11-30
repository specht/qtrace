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

#pragma once
#include <QtCore>
#include <ptb/ScanIterator.h>
#include "IsotopeEnvelope.h"


#define HYDROGEN 1.007947
#define WATER 18.01057
//#define NEUTRON 1.008664915
// TODO: verify this neutron mass
#define NEUTRON 1.002
#define HEAVY_ARGININE 6.020129027
#define HEAVY_PROLINE 5.016774189
#define HEAVY_NITROGEN 0.997034893

struct r_LabelType
{
	enum Enumeration {
		HeavyArginine = 0,
		HeavyArginineAndProline = 1,
		N15Labeling = 2,
		Size
	};
};


struct r_ScanQuantitationResult
{
	bool mb_IsGood;
    QString ms_Peptide;
	double md_AmountUnlabeled;
	double md_AmountLabeled;
	double md_Snr;
	double md_RetentionTime;
	QList<double> mk_TargetMz;
	QList<double> mk_ForbiddenMz;
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
                 bool ab_UseArea = true,
				 QList<tk_IntPair> ak_MsLevels = QList<tk_IntPair>() << tk_IntPair(0, 0x10000),
				 int ai_IsotopeCount = 3, int ai_MinCharge = 2, int ai_MaxCharge = 3, 
				 double ad_MinSnr = 2.0, double ad_MassAccuracy = 5.0,
				 double ad_ExcludeMassAccuracy = 10.0,
				 QIODevice* ak_CsvOutDevice_ = NULL, QIODevice* ak_XhtmlOutDevice_ = NULL,
				 bool ab_CheckLightForbiddenPeaks = true,
				 bool ab_CheckHeavyForbiddenPeaks = false,
				 bool ab_PrintStatusMessages = true,
                 bool ab_LogScale = true);
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
					 QList<double> ak_TargetMz, 
					 QList<double> ak_ForbiddenMz);
	void fitGaussian(double* a_, double* b_, double* c_, double x0, double y0, 
					 double x1, double y1, double x2, double y2);
	double gaussian(double x, double a, double b, double c);
    QHash<QString, int> compositionForPeptide(const QString& as_Peptide);

/*
	peptide:
		spectra file/charge:
			ratio:
			certainty:
			more info:
*/
	QTextStream mk_CsvOutStream;
	QTextStream mk_XhtmlOutStream;
	r_LabelType::Enumeration me_LabelType;
    bool mb_UseArea;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_MassAccuracy;
	double md_ExcludeMassAccuracy;
	double md_ElutionProfilePeakWidth;
	bool mb_CheckLightForbiddenPeaks;
	bool mb_CheckHeavyForbiddenPeaks;
	bool mb_PrintStatusMessages;
    bool mb_LogScale;
	QList<double> mk_AllTargetMasses;
	QString ms_CurrentSpot;
	QStringList mk_Peptides;
	QHash<char, double> mk_AminoAcidWeight;
    QHash<char, QHash<QString, int> > mk_AminoAcidComposition;
	unsigned int mui_QuantitationResultCount;
    
    QHash<QString, QSet<QString> > mk_UnlabeledRequiredTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_UnlabeledConsideredLeftTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_UnlabeledConsideredRightTargetMzForPeptideCharge;
    QHash<QString, QSet<QString> > mk_LabeledRequiredTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_LabeledConsideredLeftTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_LabeledConsideredRightTargetMzForPeptideCharge;
	
	// peptide-charge-label-isotope
	QHash<QString, int> mk_TargetMzIndex;
	
    k_IsotopeEnvelope mk_IsotopeEnvelope;
    k_IsotopeEnvelope mk_IsotopeEnvelopeHeavyN15;
    //QHash<QString, QList<double> > mk_IsotopeEnvelopesLight;
	
	//QHash<QString, QList<QList<double> > > mk_ElutionProfile;
};
