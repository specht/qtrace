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

#pragma once
#include <QtCore>
#include <ptb/ScanIterator.h>
#include <ptb/IsotopeEnvelope.h>
#include <ptb/RefPtr.h>


#define DEFAULT_LABEL "15N"
#define DEFAULT_SCAN_TYPE "all"
#define DEFAULT_USE_ISOTOPE_ENVELOPES true
#define DEFAULT_MIN_CHARGE 2
#define DEFAULT_MAX_CHARGE 3
#define DEFAULT_MIN_SNR 2.0
#define DEFAULT_MASS_ACCURACY 5.0
#define DEFAULT_REQUIRE_ABUNDANCE 0.5
#define DEFAULT_CONSIDER_ABUNDANCE 0.05
#define DEFAULT_MAX_FIT_ERROR 0.05
#define DEFAULT_CHECK_FORBIDDEN_PEAK true
#define DEFAULT_QUIET false
#define DEFAULT_CSV_OUTPUT true
#define DEFAULT_XHTML_OUTPUT true
#define DEFAULT_LOG_SCALE true

struct r_Parameter
{
    enum Enumeration
    {
        Label,
        ScanType,
        UseIsotopeEnvelopes,
        MinCharge,
        MaxCharge,
        MinSnr,
        MassAccuracy,
        RequireAbundance,
        ConsiderAbundance,
        MaxFitError,
        CheckForbiddenPeak,
        Quiet,
        CsvOutput,
        CsvOutputPath,
        XhtmlOutput,
        XhtmlOutputPath,
        LogScale
    };
};


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


struct r_EnvelopePeaks
{
    QSet<int> mk_RequiredIds;
    QSet<int> mk_ConsideredIds;
    QSet<int> mk_ForbiddenIds;
};


// tk_ArtificialEnvironment defines an artificial environment in
// which the abundances of the isotopes are modified (like 99% 15N)
typedef QHash<QString, r_IsotopeAbundance> tk_ArtificialEnvironment;

typedef QPair<int, r_Scan> tk_IntScanPair;
typedef QList<r_Scan> tk_ScanList;


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
typedef QHash<QString, double> tk_StringDoubleHash;


class k_QuantifierBase: public k_ScanIterator
{
public:
	k_QuantifierBase(QStringList& ak_Arguments, QSet<r_Parameter::Enumeration> ak_Parameters, QString as_ProgramName, QString as_AdditionalArguments = QString());
	virtual ~k_QuantifierBase();
	
	virtual void progressFunction(QString as_ScanId, bool ab_InterestingScan);
	
	virtual QString renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult);
    QHash<QString, int> compositionForPeptide(const QString& as_Peptide);
    void leastSquaresFit(QList<tk_DoublePair> ak_Pairs, double* ad_Factor_, double* ad_Error_);
    
    // parallel mass matching, attention: both lists must be sorted!
    QHash<int, int> matchTargetsToPeaks(QList<double> ak_PeakMz, QMultiMap<double, int> ak_Targets, double ad_MassAccuracy);
    QHash<int, int> extractMatches(QSet<int> ak_Ids, QList<int> ak_TargetIdsSorted, QHash<int, int> ak_Matches);
    virtual void calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_);
    
    void printUsageAndExit();
    bool stringToBool(QString as_String);
    int stringToInt(QString as_String);
    double stringToDouble(QString as_String);
    
    void removeNonPeptides(QSet<QString>& ak_List);
	
protected:
    virtual void parseArguments(QStringList& ak_Arguments);
	virtual double calculatePeptideMass(QString as_Peptide, int ai_Charge);
	virtual double scale(const double ad_Value) const;
	virtual r_ScanQuantitationResult 
		checkResult(QHash<int, r_Peak> ak_LightPeaksInclude, 
					 QHash<int, r_Peak> ak_HeavyPeaksInclude,
					 QHash<int, r_Peak> ak_LightPeaksExclude, 
					 QHash<int, r_Peak> ak_HeavyPeaksExclude,
					 r_Scan& ar_Scan, QString as_Peptide, 
					 int ai_Charge, 
					 QList<double> ak_TargetMz, 
					 QList<double> ak_ForbiddenMz);
	double gaussian(double x, double a, double b, double c);
    void parseLabel();
    QStringList tokenize(QString as_String);
    QString fetchNextToken(QStringList* ak_StringList_, QVariant::Type* ae_Type_);
    QVariant::Type peekNextToken(QStringList ak_StringList);
    tk_IsotopeEnvelope lightEnvelopeForPeptide(QString as_Peptide);
    tk_IsotopeEnvelope heavyEnvelopeForPeptide(QString as_Peptide);

    // general information
    QSet<r_Parameter::Enumeration> mk_Parameters;
    QString ms_ProgramName;
    QString ms_AdditionalArguments;

    // parameters
	QString ms_Label;
    bool mb_UseIsotopeEnvelopes;
	int mi_MinCharge;
	int mi_MaxCharge;
	double md_MinSnr;
	double md_MassAccuracy;
    double md_RequireAbundance;
    double md_ConsiderAbundance;
    double md_MaxFitError;
    bool mb_CheckForbiddenPeak;
	bool mb_Quiet;
    bool mb_LogScale;
    
	QHash<char, double> mk_AminoAcidWeight;
    QHash<char, QHash<QString, int> > mk_AminoAcidComposition;
    
    RefPtr<QIODevice> mk_pCsvDevice;
    RefPtr<QIODevice> mk_pXhtmlDevice;
    
    RefPtr<QTextStream> mk_pCsvStream;
    RefPtr<QTextStream> mk_pXhtmlStream;
    
	// peptide-charge-label-isotope
	QHash<QString, int> mk_TargetMzIndex;
	
    double md_WaterMass;
    double md_HydrogenMass;
    
    QString ms_CurrentSpectraFile;

    // natural isotope envelope generator
    k_IsotopeEnvelope mk_IsotopeEnvelope;
    // this hash contains a heavy isotope generator for every labeled amino acid
    QHash<QString, k_IsotopeEnvelope> mk_HeavyIsotopeEnvelopeForAminoAcid;
    QHash<QString, tk_ArtificialEnvironment> mk_Label;
};
