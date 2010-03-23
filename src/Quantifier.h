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
#include <ptb/IsotopeEnvelope.h>
#include "QuantifierBase.h"

typedef QPair<double, double> tk_DoublePair;


class k_Quantifier: public k_QuantifierBase
{
public:
	k_Quantifier(QStringList& ak_Arguments, QSet<r_Parameter::Enumeration> ak_Parameters, QString as_ProgramName, QString as_AdditionalArguments = QString());
	virtual ~k_Quantifier();
	
	virtual void run();
    
	virtual void handleScan(r_Scan& ar_Scan, bool& ab_Continue);
protected:
    virtual void parseArguments(QStringList& ak_Arguments);
    
    QStringList mk_SpectraFiles;
    QSet<QString> mk_Peptides;
};
