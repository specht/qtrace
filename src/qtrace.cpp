/*
Copyright (c) 2007-2008 Michael Specht

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

// define program name and parameters here

#include <QtCore>
#include <stdio.h>
#include <ptb/FastaReader.h>
#include "Quantifier.h"
#include <ptb/RefPtr.h>
#include "version.h"


int main(int ai_ArgumentCount, char** ac_Arguments__)
{
    QSet<r_Parameter::Enumeration> lk_Parameters = 
        QSet<r_Parameter::Enumeration>() 
            << r_Parameter::Label
            << r_Parameter::ScanType
            << r_Parameter::UseIsotopeEnvelopes
            << r_Parameter::MinCharge
            << r_Parameter::MaxCharge
            << r_Parameter::MinSnr
            << r_Parameter::MassAccuracy
            << r_Parameter::AbsenceMassAccuracyFactor
            << r_Parameter::RequireAbundance
            << r_Parameter::ConsiderAbundance
            << r_Parameter::MaxFitError
            << r_Parameter::CheckForbiddenPeak
            << r_Parameter::Quiet
            << r_Parameter::LogScale
            << r_Parameter::CsvOutput
            << r_Parameter::CsvOutputPath
            << r_Parameter::XhtmlOutput
            << r_Parameter::XhtmlOutputPath;

    QStringList lk_Arguments;
	for (int i = 1; i < ai_ArgumentCount; ++i)
		lk_Arguments << ac_Arguments__[i];
    
    k_Quantifier lk_Quantifier(lk_Arguments, lk_Parameters, "qtrace", "--spectraFiles <path> ... --peptideFiles <path> ... --peptides <peptide> ...");
    lk_Quantifier.run();
        
    /*
    k_Quantifier 
        lk_Quantifier(ls_Label, le_ScanType, le_AmountEstimation,
                      QList<tk_IntPair>() << tk_IntPair(1, 1),
                      li_MinCharge, li_MaxCharge, ld_MinSnr, 
                      ld_MassAccuracy, ld_RequireAbundance, 
                      ld_ConsiderAbundance, ld_MaxFitError, lk_CsvDevice_, 
                      lk_XhtmlDevice_, lb_CheckForbiddenPeak, 
                      lb_PrintStatusMessages, lb_LogScale);
                      */

                      /*
	// remove duplicate peptides
	QSet<QString> lk_PeptidesSet = lk_Peptides.toSet();
	lk_Peptides = lk_PeptidesSet.toList();
		
	lk_Quantifier.quantify(lk_SpectraFiles, lk_Peptides);
    */
}
