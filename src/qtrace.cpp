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
	printf("Usage: qtrace [options] --spectraFiles [spectra files] --peptides [peptides] --peptideFiles [peptide files]\n");
	printf("Spectra files may be mzData, mzXML or mzML, optionally compressed (.gz|.bz2|.zip).\n");
	printf("Options:\n");
	printf("  --label [R|RP|N15] (default: RP)\n");
	printf("      R  : SILAC labeling with 13C6-arginine.\n");
	printf("      RP : SILAC labeling with 13C6-arginine, additional care is taken for accidentally labeled proline residues.\n");
	printf("      N15: N15 labeling, where every amino acid is affected.\n");
	printf("  --scanType [full|sim|all] (default: all)\n");
	printf("  --isotopeCount [int] (default: 3)\n");
	printf("  --minCharge [int] (default: 2)\n");
	printf("  --maxCharge [int] (default: 3)\n");
	printf("  --minSnr [float] (default: 2.0)\n");
	printf("  --massAccuracy (ppm) [float] (default: 5.0)\n");
	printf("      This mass accuracy is used to check for the presence of peaks.\n");
	printf("  --excludeMassAccuracy (ppm) [float] (default: 30.0)\n");
	printf("      This mass accuracy is used to check for the absence of peaks.\n");
	printf("  --csvOutput [flag] (default: yes)\n");
	printf("      Enable or disable CSV output.\n");
	printf("  --csvOutputTarget [path] (default: stdout)\n");
	printf("      Redirect CSV output to a file. Enables CSV output.\n");
	printf("  --checkLightForbiddenPeaks [flag] (default: yes)\n");
	printf("      Check for light forbidden peak absence.\n");
	printf("  --checkHeavyForbiddenPeaks [flag] (default: no)\n");
	printf("      Check for heavy forbidden peak absence.\n");
/*	printf("  --xhtmlOutput [flag] (default: no)\n");
	printf("      Enable or disable XHTML output.\n");
	printf("  --xhtmlOutputTarget [string] (default: stdout)\n");
	printf("      Redirect XHTML output to a file. Enables XHTML output.\n");*/
	printf("  --statistics [flag] (default: no)\n");
	printf("      Print details about reasons why quantitation failed in scans.\n");
	printf("  --version\n");
	printf("      Print version and exit.\n");
	exit(1);
}


bool stringToBool(QString& as_String)
{
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
		printf("qtrace %s\n", gs_Version.toStdString().c_str());
		exit(0);
	}

	r_LabelType::Enumeration le_LabelType = r_LabelType::HeavyArginineAndProline;
	r_ScanType::Enumeration le_ScanType = r_ScanType::All;
	int li_IsotopeCount = 3;
	int li_MinCharge = 2;
	int li_MaxCharge = 3;
	double ld_MinSnr = 2.0;
	double ld_MassAccuracy = 5.0;
	double ld_ExcludeMassAccruracy = 30.0;
	bool lb_PrintStatistics = false;
	bool lb_CheckLightForbiddenPeaks = true;
	bool lb_CheckHeavyForbiddenPeaks = false;
	
	QFile lk_StdOut;
	lk_StdOut.open(stdout, QIODevice::WriteOnly);
	
	RefPtr<QFile> lk_pCsvOutFile;
	RefPtr<QFile> lk_pXhtmlOutFile;
	
	QIODevice* lk_CsvDevice_ = &lk_StdOut;
// 	QIODevice* lk_XhtmlDevice_ = &lk_StdOut;
	QIODevice* lk_XhtmlDevice_ = NULL;
	
	// consume options
	int li_Index;
	
	li_Index = lk_Arguments.indexOf("--label");
	if (li_Index > -1)
	{
		QString ls_Label = lk_Arguments[li_Index + 1].toUpper();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
		if (ls_Label == "R")
			le_LabelType = r_LabelType::HeavyArginine;
		else if (ls_Label == "RP")
			le_LabelType = r_LabelType::HeavyArginineAndProline;
		else if (ls_Label == "N15")
			le_LabelType = r_LabelType::N15Labeling;
		else
		{
			printf("Error: unknown label %s.\n", ls_Label.toStdString().c_str());
			exit(1);
		}
	}
	
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
	
	li_Index = lk_Arguments.indexOf("--isotopeCount");
	if (li_Index > -1)
	{
		li_IsotopeCount = QVariant(lk_Arguments[li_Index + 1]).toInt();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
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
	
	li_Index = lk_Arguments.indexOf("--excludeMassAccuracy");
	if (li_Index > -1)
	{
		ld_ExcludeMassAccruracy = QVariant(lk_Arguments[li_Index + 1]).toDouble();
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
	
	li_Index = lk_Arguments.indexOf("--checkLightForbiddenPeaks");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		lb_CheckLightForbiddenPeaks = stringToBool(ls_Value);
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--checkHeavyForbiddenPeaks");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		lb_CheckHeavyForbiddenPeaks = stringToBool(ls_Value);
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
/*	li_Index = lk_Arguments.indexOf("--xhtmlOutput");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		if (!stringToBool(ls_Value))
			lk_XhtmlDevice_ = NULL;
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
	*/
	
	li_Index = lk_Arguments.indexOf("--printStatistics");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		lb_PrintStatistics = stringToBool(ls_Value);
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	//RefPtr<QIODevice> lk_pTextDevice(new QIODevice(stdout));
	
	k_Quantifier lk_Quantifier(le_LabelType, le_ScanType,
		QList<tk_IntPair>() << tk_IntPair(1, 1),
		li_IsotopeCount, li_MinCharge, li_MaxCharge, ld_MinSnr, ld_MassAccuracy, 
		ld_ExcludeMassAccruracy, lk_CsvDevice_, lk_XhtmlDevice_, lb_PrintStatistics,
		lb_CheckLightForbiddenPeaks, lb_CheckHeavyForbiddenPeaks);
		
	QStringList lk_SpectraFiles;
	QStringList lk_Peptides;
	while (!lk_Arguments.empty())
	{
		if (lk_Arguments.first() == "--spectraFiles")
		{
			lk_Arguments.removeFirst();
			while (!lk_Arguments.empty() && !lk_Arguments.first().startsWith("-"))
				lk_SpectraFiles << lk_Arguments.takeFirst();
		}
		else if (lk_Arguments.first() == "--peptides")
		{
			lk_Arguments.removeFirst();
			while (!lk_Arguments.empty() && !lk_Arguments.first().startsWith("-"))
				lk_Peptides.push_back(lk_Arguments.takeFirst().toUpper());
		}
		else if (lk_Arguments.first() == "--peptideFiles")
		{
			lk_Arguments.removeFirst();
			while (!lk_Arguments.empty() && !lk_Arguments.first().startsWith("-"))
			{	
				QString ls_Path = lk_Arguments.takeFirst();
				QFile lk_File(ls_Path);
				if (!lk_File.open(QIODevice::ReadOnly))
				{
					printf("Error: Unable to open %s.\n", ls_Path.toStdString().c_str());
					exit(1);
				}
				QTextStream lk_Stream(&lk_File);
				while (!lk_Stream.atEnd())
				{
					QString ls_Line = lk_Stream.readLine().trimmed();
					if (!ls_Line.startsWith(">"))
						lk_Peptides.push_back(ls_Line.toUpper());
				}
			}
		}
		else
		{
			printf("Error: Unknown command line switch: %s\n", lk_Arguments.first().toStdString().c_str());
			exit(1);
		}
	}
	
	if (lk_SpectraFiles.empty() || lk_Peptides.empty())
		printUsageAndExit();
		
	// remove duplicate peptides
	QSet<QString> lk_PeptidesSet = lk_Peptides.toSet();
	lk_Peptides = lk_PeptidesSet.toList();
		
	lk_Quantifier.quantify(lk_SpectraFiles, lk_Peptides);
}
