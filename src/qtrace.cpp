/*
Copyright (c) 2007-20089 Michael Specht

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

#include <QtCore>
#include <stdio.h>
#include "Quantifier.h"
#include "RefPtr.h"
#include "version.h"


void printUsageAndExit()
{
	printf("Usage: qtrace [options] --spectraFiles [spectra files] --peptides [peptides] --peptideFiles [peptide files]\n");
	printf("Spectra files may be mzData, mzXML or mzML, optionally compressed (.gz|.bz2|.zip).\n");
	printf("Options:\n");
	printf("  --scanType [full|sim|all] (default: all)\n");
	printf("  --isotopeCount [int] (default: 3)\n");
	printf("  --minCharge [int] (default: 2)\n");
	printf("  --maxCharge [int] (default: 3)\n");
	printf("  --minSnr [float] (default: 2.0)\n");
	printf("  --massAccuracy (ppm) [float] (default: 5.0)\n");
	printf("  --svgOutPath [string] (default: none)\n");
	printf("      If svgOutPath is specified, SVG files will be created for each\n");
	printf("      specturm that was used for the quantitation. Please note that\n");
	printf("      the specified path must be empty.\n");
	printf("  --textOutput [flag] (default: yes)\n");
	printf("      Enable or disable text output.\n");
	printf("  --textOutputTarget [path] (default: stdout)\n");
	printf("      Redirect text output to a file. Enables text output.\n");
	printf("  --yamlOutput [flag] (default: no)\n");
	printf("      Enable or disable YAML output.\n");
	printf("  --yamlOutputTarget [path] (default: stdout)\n");
	printf("      Redirect YAML output to a file. Enables YAML output.\n");
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
		
	r_ScanType::Enumeration le_ScanType = r_ScanType::All;
	int li_IsotopeCount = 3;
	int li_MinCharge = 2;
	int li_MaxCharge = 3;
	double ld_MinSnr = 2.0;
	double ld_MassAccuracy = 5.0;
	bool lb_PrintStatistics = false;
	QString ls_SvgOutPath;
	
	QFile lk_StdOut;
	lk_StdOut.open(stdout, QIODevice::WriteOnly);
	
	RefPtr<QFile> lk_pTextOutFile;
	RefPtr<QFile> lk_pYamlOutFile;
	
	QIODevice* lk_TextDevice_ = &lk_StdOut;
	QIODevice* lk_YamlDevice_ = &lk_StdOut;
	
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
	
	li_Index = lk_Arguments.indexOf("--svgOutPath");
	if (li_Index > -1)
	{
		ls_SvgOutPath = lk_Arguments[li_Index + 1];
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--svgOutPath");
	if (li_Index > -1)
	{
		ls_SvgOutPath = lk_Arguments[li_Index + 1];
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--textOutput");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		if (!stringToBool(ls_Value))
			lk_TextDevice_ = NULL;
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--textOutputTarget");
	if (li_Index > -1)
	{
		lk_pTextOutFile = RefPtr<QFile>(new QFile(lk_Arguments[li_Index + 1]));
		lk_pTextOutFile->open(QIODevice::WriteOnly);
		lk_TextDevice_ = lk_pTextOutFile.get_Pointer();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--yamlOutput");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		if (!stringToBool(ls_Value))
			lk_YamlDevice_ = NULL;
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--yamlOutputTarget");
	if (li_Index > -1)
	{
		lk_pYamlOutFile = RefPtr<QFile>(new QFile(lk_Arguments[li_Index + 1]));
		lk_pYamlOutFile->open(QIODevice::WriteOnly);
		lk_YamlDevice_ = lk_pYamlOutFile.get_Pointer();
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	li_Index = lk_Arguments.indexOf("--printStatistics");
	if (li_Index > -1)
	{
		QString ls_Value = lk_Arguments[li_Index + 1];
		lb_PrintStatistics = stringToBool(ls_Value);
		lk_Arguments.removeAt(li_Index);
		lk_Arguments.removeAt(li_Index);
	}
	
	//RefPtr<QIODevice> lk_pTextDevice(new QIODevice(stdout));
	
	k_Quantifier lk_Quantifier(le_ScanType, QList<tk_IntPair>() << tk_IntPair(1, 1),
		li_IsotopeCount, li_MinCharge, li_MaxCharge, ld_MinSnr,
		ld_MassAccuracy, ls_SvgOutPath, lk_TextDevice_, lk_YamlDevice_, lb_PrintStatistics);
		
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
