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
#include <ptb/IsotopeEnvelope.h>


struct r_IsotopeAbundance
{
    r_IsotopeAbundance()
        : mi_Isotope(0)
        , mf_Efficiency(1.0)
    {
    }
    
    r_IsotopeAbundance(int ai_Isotope, float af_Efficiency)
        : mi_Isotope(ai_Isotope)
        , mf_Efficiency(af_Efficiency)
    {
    }
    
    int mi_Isotope;
    float mf_Efficiency;
};

// tk_ArtificialEnvironment defines an artificial environment in
// which the abundances of the isotopes are modified (like 99% 15N)
typedef QHash<QString, r_IsotopeAbundance> tk_ArtificialEnvironment;


struct r_AmountEstimation
{
    enum Enumeration {
        Profile = 0,
        Intensity = 1,
        Area = 2,
        Size
    };
};


struct r_ScanQuantitationResult
{
	bool mb_IsGood;
    QString ms_Peptide;
	double md_AmountUnlabeled;
	double md_AmountLabeled;
    double md_UnlabeledProfileScale;
    double md_LabeledProfileScale;
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
    double md_UnlabeledError;
    double md_LabeledError;
	
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

typedef QPair<double, double> tk_DoublePair;


class k_Quantifier: public k_ScanIterator
{
public:
	k_Quantifier(QString as_Label = "15N",
				 r_ScanType::Enumeration ae_ScanType = r_ScanType::All,
                 r_AmountEstimation::Enumeration ab_UseArea = r_AmountEstimation::Profile,
				 QList<tk_IntPair> ak_MsLevels = QList<tk_IntPair>() << tk_IntPair(0, 0x10000),
				 int ai_MinCharge = 2, int ai_MaxCharge = 3, 
				 double ad_MinSnr = 2.0, double ad_MassAccuracy = 5.0,
                 double ad_RequireAbundance = 0.5,
                 double ad_ConsiderAbundance = 0.05,
                 double ad_MaxFitError = 0.01,
				 QIODevice* ak_CsvOutDevice_ = NULL, QIODevice* ak_XhtmlOutDevice_ = NULL,
				 bool ab_CheckForbiddenPeak = true,
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
    double leastSquaresFit(QList<tk_DoublePair> ak_Pairs);
    void parseLabel();
    QStringList tokenize(QString as_String);
    QString fetchNextToken(QStringList* ak_StringList_, QVariant::Type* ae_Type_);
    QVariant::Type peekNextToken(QStringList ak_StringList);
    tk_IsotopeEnvelope heavyEnvelopeForPeptide(QString as_Peptide);

	QTextStream mk_CsvOutStream;
	QTextStream mk_XhtmlOutStream;
	QString ms_Label;
    r_AmountEstimation::Enumeration me_AmountEstimation;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_MassAccuracy;
	bool mb_CheckForbiddenPeak;
    double md_RequireAbundance;
    double md_ConsiderAbundance;
    double md_MaxFitError;
	bool mb_PrintStatusMessages;
    bool mb_LogScale;
	QList<double> mk_AllTargetMasses;
	QString ms_CurrentSpot;
	QStringList mk_Peptides;
	QHash<char, double> mk_AminoAcidWeight;
    QHash<char, QHash<QString, int> > mk_AminoAcidComposition;
    
    QHash<QString, QSet<QString> > mk_UnlabeledRequiredTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_UnlabeledConsideredLeftTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_UnlabeledConsideredRightTargetMzForPeptideCharge;
    QHash<QString, QSet<QString> > mk_LabeledRequiredTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_LabeledConsideredLeftTargetMzForPeptideCharge;
    QHash<QString, QStringList> mk_LabeledConsideredRightTargetMzForPeptideCharge;
    
    QHash<QString, tk_IsotopeEnvelope> mk_UnlabeledIsotopeEnvelopeForPeptideCharge;
    QHash<QString, tk_IsotopeEnvelope> mk_LabeledIsotopeEnvelopeForPeptideCharge;
	
	// peptide-charge-label-isotope
	QHash<QString, int> mk_TargetMzIndex;
	
    k_IsotopeEnvelope mk_IsotopeEnvelope;
    //QHash<QString, QList<double> > mk_IsotopeEnvelopesLight;
	
	//QHash<QString, QList<QList<double> > > mk_ElutionProfile;
    double md_WaterMass;
    double md_HydrogenMass;
    
    QHash<QString, k_IsotopeEnvelope> mk_HeavyIsotopeEnvelopeForAminoAcid;
    QHash<QString, tk_ArtificialEnvironment> mk_Label;
};
