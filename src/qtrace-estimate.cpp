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


void printUsageAndExit()
{
	printf("Usage: qtrace-estimate [options] [spectra file] [peptide 1] [RT 1] [peptide 2] ...\n");
	printf("The spectra file may be mzData, mzXML or mzML, optionally compressed (.gz|.bz2|.zip).\n");
	printf("Options:\n");
	printf("  --scanType [full|sim|all] (default: all)\n");
	printf("  --minCharge [int] (default: 2)\n");
	printf("  --maxCharge [int] (default: 3)\n");
	printf("  --minSnr [float] (default: 2.0)\n");
	printf("  --massAccuracy (ppm) [float] (default: 5.0)\n");
	printf("      This mass accuracy is used to check for the presence of peaks.\n");
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

	QString ls_Label = "15N";
	r_ScanType::Enumeration le_ScanType = r_ScanType::All;
    r_AmountEstimation::Enumeration le_AmountEstimation = r_AmountEstimation::Profile;
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
	
	QFile lk_StdOut;
	lk_StdOut.open(stdout, QIODevice::WriteOnly);
	
	RefPtr<QFile> lk_pCsvOutFile;
	RefPtr<QFile> lk_pXhtmlOutFile;
	
	QIODevice* lk_CsvDevice_ = &lk_StdOut;
// 	QIODevice* lk_XhtmlDevice_ = &lk_StdOut;
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
    
    lk_SpectraFiles.append(lk_Arguments.takeFirst());
    while (!lk_Arguments.empty())
    {
        if (lk_Arguments.size() >= 2)
        {
            QString ls_Peptide = lk_Arguments.takeFirst();
            QString ls_Rt = lk_Arguments.takeFirst();
            lk_Peptides[ls_Peptide] = ls_Rt.toDouble();
        }
    }
    
	if (lk_SpectraFiles.empty() || lk_Peptides.empty())
		printUsageAndExit();
		
    k_Quantifier 
        lk_Quantifier(ls_Label, le_ScanType, le_AmountEstimation,
                      QList<tk_IntPair>() << tk_IntPair(1, 1),
                      li_MinCharge, li_MaxCharge, ld_MinSnr, 
                      ld_MassAccuracy, ld_RequireAbundance, 
                      ld_ConsiderAbundance, ld_MaxFitError, lk_CsvDevice_, 
                      lk_XhtmlDevice_, lb_CheckForbiddenPeak, 
                      lb_PrintStatusMessages, lb_LogScale);
        
	lk_Quantifier.quantify(lk_SpectraFiles, lk_Peptides.keys());
}
