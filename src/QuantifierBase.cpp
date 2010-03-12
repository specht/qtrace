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

#include "QuantifierBase.h"
#include <QtCore>
#include <QtSvg>
#include <math.h> 
#include <limits>
#include "Tango.h"

#define ABSENCE_MASS_ACCURACY_FACTOR 2.0


extern QString gs_Version;


k_QuantifierBase::k_QuantifierBase(QStringList& ak_Arguments, QSet<r_Parameter::Enumeration> ak_Parameters, QString as_ProgramName, QString as_AdditionalArguments)
    : k_ScanIterator()
    , mk_Parameters(ak_Parameters)
    , ms_ProgramName(as_ProgramName)
    , ms_AdditionalArguments(as_AdditionalArguments)
	, ms_Label(DEFAULT_LABEL)
    , mb_UseIsotopeEnvelopes(DEFAULT_USE_ISOTOPE_ENVELOPES)
	, mi_MinCharge(DEFAULT_MIN_CHARGE)
	, mi_MaxCharge(DEFAULT_MAX_CHARGE)
	, md_MinSnr(DEFAULT_MIN_SNR)
	, md_MassAccuracy(DEFAULT_MASS_ACCURACY)
    , md_RequireAbundance(DEFAULT_REQUIRE_ABUNDANCE)
    , md_ConsiderAbundance(DEFAULT_CONSIDER_ABUNDANCE)
    , md_MaxFitError(DEFAULT_MAX_FIT_ERROR)
	, mb_CheckForbiddenPeak(DEFAULT_CHECK_FORBIDDEN_PEAK)
	, mb_Quiet(DEFAULT_QUIET)
    , mb_LogScale(DEFAULT_LOG_SCALE)
{
    Q_INIT_RESOURCE(qtrace);
    
    setMsLevels(QList<tk_IntPair>() << tk_IntPair(1, 1));
    parseArguments(ak_Arguments);
    
    md_HydrogenMass = mk_IsotopeEnvelope.mk_BaseIsotopeMass["H"];
    double ld_OxygenMass = mk_IsotopeEnvelope.mk_BaseIsotopeMass["O"];
    md_WaterMass = md_HydrogenMass * 2.0 + ld_OxygenMass;
    
	// test parameter sanity
	if (mi_MinCharge < 1)
	{
		printf("Error: minimum charge (%d) out of range.\n", mi_MinCharge);
		exit(1);
	}
	if (mi_MaxCharge < 1)
	{
		printf("Error: maximum charge (%d) out of range.\n", mi_MaxCharge);
		exit(1);
	}
	if (mi_MaxCharge < mi_MinCharge)
	{
		printf("Error: minimum charge (%d) bigger than maximum charge(%d).\n", mi_MinCharge, mi_MaxCharge);
		exit(1);
	}
	if (md_MassAccuracy < 0.0)
	{
		printf("Error: mass accuracy must not be less than zero.\n");
		exit(1);
	}
    if (md_RequireAbundance < 0.0 || md_RequireAbundance > 1.0)
    {
        printf("Error: require abundance out of range (0.0 - 1.0).\n");
        exit(1);
    }
    if (md_ConsiderAbundance < 0.0 || md_ConsiderAbundance > 1.0)
    {
        printf("Error: consider abundance out of range (0.0 - 1.0).\n");
        exit(1);
    }
    if (md_ConsiderAbundance > md_RequireAbundance)
    {
        printf("Error: consider abundance must not be greater than require abundance.\n");
        exit(1);
    }
	
	QFile lk_File(":ext/proteomics-knowledge-base/amino-acids.csv");
	if (!lk_File.open(QFile::ReadOnly))
    {
        fprintf(stderr, "Error: Unable to open amino-acids.csv.\n");
        exit(1);
    }
	QTextStream lk_TextStream(&lk_File);

    QString ls_Header = lk_TextStream.readLine().toLower().trimmed();
    QStringList lk_Header = ls_Header.split(",");
    int li_LetterIndex = lk_Header.indexOf("single letter code");
    int li_CompositionIndex = lk_Header.indexOf("composition");
    int li_MassIndex = lk_Header.indexOf("monoisotopic mass");
    if (li_LetterIndex < 0 || li_CompositionIndex < 0 || li_MassIndex < 0)
    {
        fprintf(stderr, "Error: Something is wrong with amino-acids.csv.\n");
        exit(1);
    }
	while (!lk_TextStream.atEnd())
	{
		QString ls_Line = lk_TextStream.readLine().trimmed();
		if (ls_Line.trimmed().isEmpty())
			continue;

		QStringList lk_List;
		foreach (QString ls_Entry, ls_Line.split(QChar(',')))
		{
            if (ls_Entry.length() >= 2)
            {
                if (ls_Entry.at(0) == QChar('"') && ls_Entry.at(ls_Entry.length() - 1) == QChar('"'))
                    ls_Entry = ls_Entry.mid(1, ls_Entry.length() - 2);
            }
			lk_List << ls_Entry;
		}
		char lc_AminoAcid = lk_List[li_LetterIndex][0].toAscii();
		mk_AminoAcidWeight[lc_AminoAcid] = lk_List[li_MassIndex].toDouble();
		QString ls_Composition = lk_List[li_CompositionIndex];
		ls_Composition.remove("\"");
		QHash<QString, int> lk_ElementCount;
		while (!ls_Composition.isEmpty())
		{
			char lc_Element = ls_Composition.at(0).toAscii();
			lk_ElementCount[QString(lc_Element)] = 1;
			ls_Composition = ls_Composition.right(ls_Composition.length() - 1);
			if (ls_Composition.isEmpty())
				break;
			QString ls_Number;
			char lc_NextChar = ls_Composition.at(0).toAscii();
			while (lc_NextChar >= '0' && lc_NextChar <= '9')
			{
				ls_Number += QChar(lc_NextChar);
				ls_Composition = ls_Composition.right(ls_Composition.length() - 1);
				if (ls_Composition.isEmpty())
					break;
				lc_NextChar = ls_Composition.at(0).toAscii();
			}
			if (!ls_Number.isEmpty())
				lk_ElementCount[QString(lc_Element)] = ls_Number.toInt();
		}
        mk_AminoAcidComposition[lc_AminoAcid] = lk_ElementCount;
	}
	lk_File.close();
    
    if (mk_pCsvDevice)
        mk_pCsvStream = RefPtr<QTextStream>(new QTextStream(mk_pCsvDevice.get_Pointer()));
    if (mk_pXhtmlDevice)
        mk_pXhtmlStream = RefPtr<QTextStream>(new QTextStream(mk_pXhtmlDevice.get_Pointer()));
    
    parseLabel();
}


k_QuantifierBase::~k_QuantifierBase()
{
}


bool compareScanQuantitationBySnr(const r_ScanQuantitationResult& a, const r_ScanQuantitationResult& b)
{
	return a.md_Snr > b.md_Snr;
}


bool compareScanQuantitationByRetentionTime(const r_ScanQuantitationResult& a, const r_ScanQuantitationResult& b)
{
	return a.md_RetentionTime < b.md_RetentionTime;
}


bool sortByMz(const QPair<double, QString>& a, const QPair<double, QString>& b)
{
	return a.first < b.first;
}


QString dtos(double ad_Value)
{
	static char lc_String_[2048];
	sprintf(lc_String_, "%1.2f", ad_Value);
	return QString(lc_String_);
}


void k_QuantifierBase::progressFunction(QString as_ScanId, bool)
{
	if (!mb_Quiet)
		printf("\r%s: scan #%s...", ms_CurrentSpectraFile.toStdString().c_str(), as_ScanId.toStdString().c_str());
}


void k_QuantifierBase::parseArguments(QStringList& ak_Arguments)
{
    if (!ak_Arguments.empty() && ak_Arguments.first() == "--version")
    {
        printf("%s %s\n", ms_ProgramName.toStdString().c_str(), gs_Version.toStdString().c_str());
        exit(0);
    }
    
    if (!ak_Arguments.empty() && ak_Arguments.first() == "--help")
        printUsageAndExit();
    
    QFile* lk_StdOut_ = new QFile();
    lk_StdOut_->open(stdout, QIODevice::WriteOnly);
    mk_pCsvDevice = RefPtr<QIODevice>(lk_StdOut_);

    // consume options
    while (!ak_Arguments.empty())
    {
        QString ls_Key = ak_Arguments.takeFirst();
        if (mk_Parameters.contains(r_Parameter::Label) && ls_Key == "--label")
            ms_Label = ak_Arguments.takeFirst();
        else if (mk_Parameters.contains(r_Parameter::ScanType) && ls_Key == "--scanType")
        {
            QString ls_Value = ak_Arguments.takeFirst();
            if (ls_Value == "sim")
                me_ScanType = r_ScanType::SIM;
            else if (ls_Value == "full")
            // why would we say MS1 and MSn here? Because I saw MS1 full scans that
            // were marked as MSn with ms level 1, which probably should be MS1 in the
            // first place, but what can a girl do? It really just means don't use SIM.
                me_ScanType = (r_ScanType::Enumeration)(r_ScanType::MS1 | r_ScanType::MSn);
            else if (ls_Value == "all")
                me_ScanType = r_ScanType::All;
            else
            {
                printf("Error: unknown scan type %s.\n", ls_Value.toStdString().c_str());
                exit(1);
            }
        }
        else if (mk_Parameters.contains(r_Parameter::CsvOutput) && ls_Key == "--csvOutput")
        {
            if (!stringToBool(ak_Arguments.takeFirst()))
                mk_pCsvDevice = RefPtr<QIODevice>();
        }
        else if (mk_Parameters.contains(r_Parameter::CsvOutputPath) && ls_Key == "--csvOutputPath")
        {
            QString ls_Path = ak_Arguments.takeFirst();
            QFile* lk_File_ = new QFile(ls_Path);
            if (!lk_File_->open(QIODevice::WriteOnly))
            {
                printf("Error: Unable to open %s for writing.\n", ls_Path.toStdString().c_str());
                exit(1);
            }
            mk_pCsvDevice = RefPtr<QIODevice>(lk_File_);
        }
        else if (mk_Parameters.contains(r_Parameter::XhtmlOutputPath) && ls_Key == "--xhtmlOutputPath")
        {
            QString ls_Path = ak_Arguments.takeFirst();
            QFile* lk_File_ = new QFile(ls_Path);
            if (!lk_File_->open(QIODevice::WriteOnly))
            {
                printf("Error: Unable to open %s for writing.\n", ls_Path.toStdString().c_str());
                exit(1);
            }
            mk_pXhtmlDevice = RefPtr<QIODevice>(lk_File_);
        }
        else if (mk_Parameters.contains(r_Parameter::UseIsotopeEnvelopes) && ls_Key == "--useIsotopeEnvelopes")
            mb_UseIsotopeEnvelopes = stringToBool(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MinCharge) && ls_Key == "--minCharge")
            mi_MinCharge = stringToInt(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MaxCharge) && ls_Key == "--maxCharge")
            mi_MaxCharge = stringToInt(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MinSnr) && ls_Key == "--minSnr")
            md_MinSnr = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MassAccuracy) && ls_Key == "--massAccuracy")
            md_MassAccuracy = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::RequireAbundance) && ls_Key == "--requireAbundance")
            md_RequireAbundance = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::ConsiderAbundance) && ls_Key == "--considerAbundance")
            md_ConsiderAbundance = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MaxFitError) && ls_Key == "--maxFitError")
            md_MaxFitError = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::CheckForbiddenPeak) && ls_Key == "--checkForbiddenPeak")
            mb_CheckForbiddenPeak = stringToBool(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::Quiet) && ls_Key == "--quiet")
            mb_Quiet = stringToBool(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::LogScale) && ls_Key == "--logScale")
            mb_LogScale = stringToBool(ak_Arguments.takeFirst());
        else
        {
            // this parameter is unknown, push it back and stop consuming options
            ak_Arguments.push_front(ls_Key);
            break;
        }
    }

/*
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
    */
}


double k_QuantifierBase::calculatePeptideMass(QString as_Peptide, int ai_Charge)
{
	double ld_Mass = md_WaterMass;
	for (int i = 0; i < as_Peptide.length(); ++i)
		ld_Mass += mk_AminoAcidWeight[as_Peptide.at(i).toAscii()];
	// UNSURE ABOUT THIS VALUE BUT OK WITH BIANCAS EXAMPLES, should be 1.0078250
	return (ld_Mass + md_HydrogenMass * ai_Charge) / ai_Charge;
}


QString k_QuantifierBase::renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult)
{
    /*
	double ld_Ratio = 4.0;
	double ld_Width = 950.0;
	double ld_Height = ld_Width / ld_Ratio;
	double ld_BorderTop = 4.0;
	double ld_BorderRight = 16.0;
	double ld_BorderBottom = 4.0;
	double ld_BorderLeft = 16.0;
	
	if (ar_QuantitationResult.md_MinMz == 0.0)
		ar_QuantitationResult.md_MinMz = ar_Scan.mr_Spectrum.md_MzValues_[0];
	if (ar_QuantitationResult.md_MaxMz == 0.0)
		ar_QuantitationResult.md_MaxMz = ar_Scan.mr_Spectrum.md_MzValues_[ar_Scan.mr_Spectrum.mi_PeaksCount - 1];
		
	int li_Start = 0;
	int li_End = ar_Scan.mr_Spectrum.mi_PeaksCount - 1;
	
	// adjust start and end values
	while (li_Start < ar_Scan.mr_Spectrum.mi_PeaksCount - 1 && ar_Scan.mr_Spectrum.md_MzValues_[li_Start] < ar_QuantitationResult.md_MinMz)
		++li_Start;
	while (li_End > 0 && ar_Scan.mr_Spectrum.md_MzValues_[li_End] > ar_QuantitationResult.md_MaxMz)
		--li_End;
	
	if (li_Start > 0)
		--li_Start;
	if (li_End < ar_Scan.mr_Spectrum.mi_PeaksCount - 1)
		++li_End;
	if (li_End < ar_Scan.mr_Spectrum.mi_PeaksCount - 1)
		++li_End;
		
	if (li_Start > li_End)
	{
		// oops! start and end got crossed.
		li_Start = 0;
		li_End = 0;
	}
	
	// determine x and y extensions
	double xmin = ar_Scan.mr_Spectrum.md_MzValues_[li_Start];
	double xmax = ar_Scan.mr_Spectrum.md_MzValues_[li_End];
	xmin = ar_QuantitationResult.md_MinMz;
	xmax = ar_QuantitationResult.md_MaxMz;
	double ymax = ar_Scan.mr_Spectrum.md_IntensityValues_[li_Start];
	for (int i = li_Start + 1; i <= li_End; ++i)
		ymax = std::max<double>(ymax, ar_Scan.mr_Spectrum.md_IntensityValues_[i]);
	// ymax remains unscaled!
	double ymaxScaled = this->scale(1.0);
    
    double ymaxUnlabeled = 0.0;
    foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_UnlabeledPeaks)
        ymaxUnlabeled = std::max<double>(ymaxUnlabeled, lr_Peak.md_PeakIntensity / ymax);
	
	QBuffer lk_Buffer;
	lk_Buffer.open(QBuffer::WriteOnly);

    QList<double> lk_Lines;
    if (mb_LogScale)
        lk_Lines << 0.01 << 0.1 << 1.0;
    else
        lk_Lines << 0.2 << 0.4 << 0.6 << 0.8 << 1.0;
            
	double x0 = ld_BorderLeft;
	double y0 = ld_Height - ld_BorderBottom;
	double dx = (ld_Width - ld_BorderLeft - ld_BorderRight) / (xmax - xmin);
	double dy = -(ld_Height - ld_BorderTop - ld_BorderBottom);
    QString ls_Labels;
	
	{
		QSvgGenerator lk_Generator;
		lk_Generator.setSize(QSize((int)ld_Width, (int)ld_Height));
		lk_Generator.setOutputDevice(&lk_Buffer);
		
		QPainter lk_Painter(&lk_Generator);
		
		// fill background
		lk_Painter.fillRect(QRectF(0.0, 0.0, ld_Width, ld_Height), QBrush(Qt::white));
		
		// draw frame
		QPen lk_Pen;
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::BevelJoin);
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Painter.setPen(lk_Pen);
		
		lk_Painter.drawLine(QPointF(0.0, y0), QPointF(ld_Width, y0));
		//lk_Painter.drawLine(QPointF(0.0, y0), QPointF(ld_Width, y0));
		//lk_Painter.drawLine(QPointF(x0, y0), QPointF(x0, y0 + (ymax - ymin) * dy));
        
        QString ls_Peptide = ar_QuantitationResult.ms_Peptide;
        int li_Charge = ar_QuantitationResult.mi_Charge;
        QString ls_PeptideChargeKey = QString("%1-%2").arg(ls_Peptide).arg(li_Charge);

        tk_IsotopeEnvelope lk_IsotopeEnvelope;
        bool lb_DrawnOnePoint;
        double ld_LastMz;
        
        lk_IsotopeEnvelope = mk_UnlabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey];
        double ld_PeptideChargeBaseMz = calculatePeptideMass(ls_Peptide, li_Charge);
        lb_DrawnOnePoint = false;
        ld_LastMz = 0.0;
        QPainterPath lk_Path;
        for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
        {
            double ld_Abundance = lk_IsotopeEnvelope[i].first * ar_QuantitationResult.md_UnlabeledProfileScale;
            double ld_Mz = lk_IsotopeEnvelope[i].second / li_Charge;
            ld_Mz += ld_PeptideChargeBaseMz;
            double x = (ld_Mz - xmin) * dx + x0;
            double y = this->scale(ld_Abundance / ymax) / ymaxScaled;
            if (y > 0.001)
            {
                y = y * dy + y0;
                if (!lb_DrawnOnePoint)
                {
                    lk_Path.moveTo(QPointF(((ld_Mz - 0.5 / li_Charge) - xmin) * dx + x0, y0));
                    lb_DrawnOnePoint = true;
                    //lk_Path.quadTo(QPointF(((ld_Mz / li_Charge) - xmin) * dx + x0, y0), QPointF(x, y));
                }
                else
                {
                    //lk_Path.lineTo(QPointF(x, y));
                }
                lk_Path.lineTo(QPointF(x, y));
                ld_LastMz = ld_Mz;
            }
        }
        lk_Path.lineTo(QPointF(((ld_LastMz + 0.5 / li_Charge) - xmin) * dx + x0, y0));
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::BevelJoin);
        lk_Pen.setColor(QColor(128, 128, 128));
        lk_Painter.setPen(lk_Pen);
        lk_Painter.setBrush(QBrush(QColor(224, 224, 224, 128)));
        
        lk_Painter.drawPath(lk_Path);
		
        lk_IsotopeEnvelope = mk_LabeledIsotopeEnvelopeForPeptideCharge[ls_PeptideChargeKey];
        lb_DrawnOnePoint = false;
        ld_LastMz = 0.0;
        lk_Path = QPainterPath();
        for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
        {
            double ld_Abundance = lk_IsotopeEnvelope[i].first * ar_QuantitationResult.md_LabeledProfileScale;
            double ld_Mz = lk_IsotopeEnvelope[i].second / li_Charge;
            ld_Mz += ld_PeptideChargeBaseMz;
            double x = (ld_Mz - xmin) * dx + x0;
            double y = this->scale(ld_Abundance / ymax) / ymaxScaled;
            if (y > 0.001)
            {
                y = y * dy + y0;
                if (!lb_DrawnOnePoint)
                {
                    lk_Path.moveTo(QPointF(((ld_Mz - 0.5 / li_Charge) - xmin) * dx + x0, y0));
                    lb_DrawnOnePoint = true;
                    //lk_Path.quadTo(QPointF(((ld_Mz / li_Charge) - xmin) * dx + x0, y0), QPointF(x, y));
                }
                else
                {
                    //lk_Path.lineTo(QPointF(x, y));
                }
                lk_Path.lineTo(QPointF(x, y));
                ld_LastMz = ld_Mz;
            }
        }
        lk_Path.lineTo(QPointF(((ld_LastMz + 0.5 / li_Charge) - xmin) * dx + x0, y0));
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::BevelJoin);
        lk_Pen.setColor(QColor(128, 128, 128));
        lk_Painter.setPen(lk_Pen);
        lk_Painter.setBrush(QBrush(QColor(224, 224, 224, 128)));
        
        lk_Painter.drawPath(lk_Path);
        
        lk_Painter.setBrush(Qt::NoBrush);
		// draw spectrum
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::RoundJoin);
		lk_Pen.setColor(QColor(TANGO_CHAMELEON_2));
		lk_Pen.setStyle(Qt::DashLine);
		lk_Painter.setPen(lk_Pen);
		
		// draw target m/z lines
		foreach (double ld_Mz, ar_QuantitationResult.mk_TargetMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		
		lk_Pen.setColor(QColor(TANGO_SCARLET_RED_1));
		lk_Painter.setPen(lk_Pen);
		// draw forbidden m/z lines
		foreach (double ld_Mz, ar_QuantitationResult.mk_ForbiddenMz)
		{
			double x = (ld_Mz - xmin) * dx + x0;
			lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
		}
		
		lk_Pen.setColor(QColor(192, 192, 192));
		lk_Painter.setPen(lk_Pen);
		
        foreach (double ld_Line, lk_Lines)
        {
            double y;
            y = this->scale(ld_Line) / ymaxScaled;
            y = y * dy + y0;
            lk_Painter.drawLine(QPointF(0.0, y), QPointF(ld_Width, y));
            ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana;' transform='rotate(90 %1 %2)'>%3%</text>\n").arg(0.0).arg(y + 3.0).arg((int)round(ld_Line * 100.0));
        }
		
		lk_Pen.setWidthF(1.0);
		lk_Pen.setJoinStyle(Qt::RoundJoin);
		lk_Pen.setColor(QColor(128, 128, 128));
		lk_Pen.setStyle(Qt::SolidLine);
		lk_Painter.setPen(lk_Pen);
		
		QVector<QPointF> lk_Points;
		for (int i = li_Start; i <= li_End; ++i)
			lk_Points.append(QPointF((ar_Scan.mr_Spectrum.md_MzValues_[i] - xmin) * dx + x0, 
									this->scale(ar_Scan.mr_Spectrum.md_IntensityValues_[i] / ymax) / ymaxScaled * dy + y0));
		lk_Painter.drawPolyline(QPolygonF(lk_Points));
		
		foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_UnlabeledPeaks)
		{
			// draw Gaussian
			lk_Pen.setWidthF(1.0);
			lk_Pen.setJoinStyle(Qt::RoundJoin);
			lk_Pen.setColor(QColor(ar_QuantitationResult.md_UnlabeledError <= md_MaxFitError ? TANGO_SKY_BLUE_1 : TANGO_SCARLET_RED_1));
			lk_Pen.setStyle(Qt::SolidLine);
			lk_Painter.setPen(lk_Pen);
			
			QVector<QPointF> lk_Points;
			double w = sqrt(lr_Peak.md_GaussC) / 2.0;
			for (double x = lr_Peak.md_PeakMz - w; x <= lr_Peak.md_PeakMz + w; x += w * 0.05)
				lk_Points.append(QPointF((x - xmin) * dx + x0, 
										this->scale(gaussian(x, lr_Peak.md_GaussA, lr_Peak.md_GaussB, lr_Peak.md_GaussC) / ymax) / ymaxScaled * dy + y0));
			lk_Painter.drawPolyline(QPolygonF(lk_Points));
		}
        
        foreach (r_Peak lr_Peak, ar_QuantitationResult.mk_LabeledPeaks)
        {
            // draw Gaussian
            lk_Pen.setWidthF(1.0);
            lk_Pen.setJoinStyle(Qt::RoundJoin);
            lk_Pen.setColor(QColor(ar_QuantitationResult.md_LabeledError <= md_MaxFitError ? TANGO_SKY_BLUE_1 : TANGO_SCARLET_RED_1));
            lk_Pen.setStyle(Qt::SolidLine);
            lk_Painter.setPen(lk_Pen);
            
            QVector<QPointF> lk_Points;
            double w = sqrt(lr_Peak.md_GaussC) / 2.0;
            for (double x = lr_Peak.md_PeakMz - w; x <= lr_Peak.md_PeakMz + w; x += w * 0.05)
                lk_Points.append(QPointF((x - xmin) * dx + x0, 
                                        this->scale(gaussian(x, lr_Peak.md_GaussA, lr_Peak.md_GaussB, lr_Peak.md_GaussC) / ymax) / ymaxScaled * dy + y0));
            lk_Painter.drawPolyline(QPolygonF(lk_Points));
        }
	}
	lk_Buffer.close();
	lk_Buffer.open(QBuffer::ReadOnly);
	lk_Buffer.seek(0);
	QString ls_Result = QString(lk_Buffer.readAll());
	lk_Buffer.close();
	foreach (r_Peak lr_Peak, (ar_QuantitationResult.mk_UnlabeledPeaks + ar_QuantitationResult.mk_LabeledPeaks))
	{
		double ld_X, ld_Y;
		ld_X = (lr_Peak.md_PeakMz - xmin) * dx + x0;
		ld_Y = this->scale(lr_Peak.md_PeakIntensity / ymax) / ymaxScaled * dy + y0;
		char lc_Number_[1024];
        lc_Number_[0] = 0;
//         if (me_AmountEstimation == r_AmountEstimation::Area)
//             sprintf(lc_Number_, "%.2g", lr_Peak.md_PeakArea);
//         else if (me_AmountEstimation == r_AmountEstimation::Intensity)
        // :TODO: now this always uses peak height, make this an option (area/intensity)
        sprintf(lc_Number_, "%.2g", lr_Peak.md_PeakIntensity);
        if (strlen(lc_Number_) > 0)
            ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana; fill: #aaa;' transform='rotate(90 %1 %2)'>%3</text>\n").arg(ld_X + 3.0).arg(ld_Y + 3.0).arg(lc_Number_);
	}
	ls_Result.replace("</svg>", ls_Labels + "</svg>");
	return ls_Result;
    */
    return QString();
}


double k_QuantifierBase::scale(const double ad_Value) const
{
    if (mb_LogScale)
        return log(ad_Value * 100.0 + 1.0);
    else
        return ad_Value;
}


void k_QuantifierBase::calculateMeanAndStandardDeviation(QList<double> ak_Values, double* ad_Mean_, double* ad_StandardDeviation_)
{
	double ld_Mean = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_Mean += ld_Value;
		
	ld_Mean /= ak_Values.size();
	
	double ld_StandardDeviation = 0.0;
	foreach (double ld_Value, ak_Values)
		ld_StandardDeviation += pow(ld_Value - ld_Mean, 2.0);
		
	ld_StandardDeviation /= ak_Values.size();
	ld_StandardDeviation = sqrt(ld_StandardDeviation);
	
	*ad_Mean_ = ld_Mean;
	*ad_StandardDeviation_ = ld_StandardDeviation;
}


void k_QuantifierBase::printUsageAndExit()
{
    //printf("Usage: %s [options] --spectraFiles [spectra files] --peptides [peptides] --peptideFiles [peptide files]\n");
    printf("Usage: %s [options] %s\n", ms_ProgramName.toStdString().c_str(), ms_AdditionalArguments.toStdString().c_str());
    printf("Spectra files may be mzData, mzXML or mzML, optionally compressed (.gz|.bz2|.zip).\n");
    printf("Options:\n");
    if (mk_Parameters.contains(r_Parameter::Label))
    {
        printf("  --label <string> (default: %s)\n", DEFAULT_LABEL);
        printf("      Specify the label here. Examples:\n");
        printf("      15N: 15N in all amino acids\n");
        printf("      R 13C: 13C in all arginine residues\n");
        printf("      RP 13C: 13C in all arginine and proline residues\n");
        printf("      RPK 13C 15N: 13C and 15N in all arginine, proline and lysine residues\n");
        printf("      R 13C K 15N: 13C in arginine, 15N in lysine residues\n");
        printf("      15N (0.99): 99%% 15N in all amino acids (remaining 1%% is N14)\n");
        printf("      The exact syntax is:\n");
        printf("        ((amino acid)* (isotope)+ (labeling efficiency)? )+\n");
        printf("      Notes:\n");
        printf("        Labeling efficiencies can only be specified for C and N.\n");
        printf("        X denotes any amino acid: 'X 15N' is the same as '15N'.\n");
    }
    if (mk_Parameters.contains(r_Parameter::ScanType))
    {
        printf("  --scanType [full|sim|all] (default: all)\n");
    }
    if (mk_Parameters.contains(r_Parameter::UseIsotopeEnvelopes))
    {
        printf("  --useIsotopeEnvelopes [yes|no] (default: %s)\n", DEFAULT_USE_ISOTOPE_ENVELOPES ? "yes" : "no");
        printf("      Specify whether amounts should be estimated via isotope envelope fitting\n");
        printf("      or whether a constant number of peaks should be used instead.\n");
    }
    if (mk_Parameters.contains(r_Parameter::MinCharge))
    {
        printf("  --minCharge <int> (default: %d)\n", DEFAULT_MIN_CHARGE);
    }
    if (mk_Parameters.contains(r_Parameter::MaxCharge))
    {
        printf("  --maxCharge <int> (default: %d)\n", DEFAULT_MAX_CHARGE);
    }
    if (mk_Parameters.contains(r_Parameter::MinSnr))
    {
        printf("  --minSnr <float> (default: %1.1f)\n", DEFAULT_MIN_SNR);
    }
    if (mk_Parameters.contains(r_Parameter::MassAccuracy))
    {
        printf("  --massAccuracy (ppm) <float> (default: %1.1f)\n", DEFAULT_MASS_ACCURACY);
        printf("      This mass accuracy is used to check for the presence of peaks.\n");
    }
    if (mk_Parameters.contains(r_Parameter::RequireAbundance))
    {
        printf("  --requireAbundance <float> (default: %1.1f)\n", DEFAULT_REQUIRE_ABUNDANCE);
        printf("      Specify which peaks must be present in an isotope envelope.\n");
    }
    if (mk_Parameters.contains(r_Parameter::ConsiderAbundance))
    {
        printf("  --considerAbundance <float> (default: %1.2f)\n", DEFAULT_CONSIDER_ABUNDANCE);
        printf("      Specify which peaks may additionally be taken into account.\n");
    }
    if (mk_Parameters.contains(r_Parameter::MaxFitError))
    {
        printf("  --maxFitError <float> (default: %1.2f)\n", DEFAULT_MAX_FIT_ERROR);
        printf("      Specify the maximum error allowed when fitting isotope envelopes.\n");
    }
    if (mk_Parameters.contains(r_Parameter::CheckForbiddenPeak))
    {
        printf("  --checkForbiddenPeak [yes|no] (default: %s)\n", DEFAULT_CHECK_FORBIDDEN_PEAK ? "yes" : "no");
        printf("      Specify whether the forbidden peak is required to be absent.\n");
    }
    if (mk_Parameters.contains(r_Parameter::CsvOutput))
    {
        printf("  --csvOutput [yes|no] (default: yes)\n");
        printf("      Enable or disable CSV output.\n");
    }
    if (mk_Parameters.contains(r_Parameter::CsvOutputPath))
    {
        printf("  --csvOutputPath <path> (default: stdout)\n");
        printf("      Redirect CSV output to a file. Enables CSV output.\n");
    }
    if (mk_Parameters.contains(r_Parameter::XhtmlOutputPath))
    {
        printf("  --xhtmlOutputPath <path> (default: none)\n");
        printf("      Output quantitation events as XHTML-embedded SVG to [path].\n");
    }
    if (mk_Parameters.contains(r_Parameter::LogScale))
    {
        printf("  --logScale [yes|no] (default: %s)\n", DEFAULT_LOG_SCALE ? "yes" : "no");
        printf("      Use logarithmic scale in XHTML spectra.\n");
    }
    if (mk_Parameters.contains(r_Parameter::Quiet))
    {
        printf("  --quiet [yes|no] (default: %s)\n", DEFAULT_QUIET ? "yes" : "no");
        printf("      Don't print status messages.\n");
    }
    printf("  --version\n");
    printf("      Print version and exit.\n");
    exit(1);
}


bool k_QuantifierBase::stringToBool(QString as_String)
{
    as_String = as_String.toLower().trimmed();
    if (as_String == "yes")
        return true;
    else if (as_String == "no")
        return false;
    else
    {
        printf("Error: Invalid flag: %s.\n", as_String.toStdString().c_str());
        exit(1);
    }
}


int k_QuantifierBase::stringToInt(QString as_String)
{
    bool lb_Ok = true;
    int li_Result=  as_String.toInt(&lb_Ok);
    if (!lb_Ok)
    {
        printf("Error: Invalid integer: %s.\n", as_String.toStdString().c_str());
        exit(1);
    }
    return li_Result;
}


double k_QuantifierBase::stringToDouble(QString as_String)
{
    bool lb_Ok = true;
    double ld_Result =  as_String.toDouble(&lb_Ok);
    if (!lb_Ok)
    {
        printf("Error: Invalid float: %s.\n", as_String.toStdString().c_str());
        exit(1);
    }
    return ld_Result;
}


void k_QuantifierBase::removeNonPeptides(QSet<QString>& ak_List)
{
    QList<QString> lk_RemovedPeptides;
    int li_RemovedPeptides = 0;
    
    foreach (QString ls_Peptide, ak_List)
    {
        QString ls_Rest = ls_Peptide;
        ls_Rest.replace(QRegExp("[GASPVTCLINDQKEMHFRYW]"), "");
        if (!ls_Rest.isEmpty())
        {
            ak_List.remove(ls_Peptide);
            ++li_RemovedPeptides;
            if (lk_RemovedPeptides.size() < 10)
                lk_RemovedPeptides.append(ls_Peptide);
        }
    }
    
    if ((!lk_RemovedPeptides.empty()) && (!mb_Quiet))
    {
        printf("Warning: Removed %d invalid peptides:\n", li_RemovedPeptides);
        foreach (QString ls_Peptide, lk_RemovedPeptides)
            printf("%s\n", ls_Peptide.toStdString().c_str());
        if (li_RemovedPeptides > lk_RemovedPeptides.size())
            printf("(more invalid peptides omitted)\n");
    }
}


r_ScanQuantitationResult 
k_QuantifierBase::checkResult(QHash<int, r_Peak> ak_LightPeaksInclude, 
						   QHash<int, r_Peak> ak_HeavyPeaksInclude,
						   QHash<int, r_Peak> ak_LightPeaksExclude, 
						   QHash<int, r_Peak> ak_HeavyPeaksExclude,
						   r_Scan& ar_Scan, QString as_Peptide, int ai_Charge,
						   QList<double> ak_TargetMz,
						   QList<double> ak_ForbiddenMz)
{
	r_ScanQuantitationResult lr_Result;
	lr_Result.mb_IsGood = false;
	
	// check whether all inclusion peaks have a good SNR
	QList<r_Peak> lk_Peaks = ak_LightPeaksInclude.values() + ak_HeavyPeaksInclude.values();
	
	// return if no inclusion peaks have been found
	if (lk_Peaks.empty())
		return lr_Result;
	
	double ld_MinSnr = lk_Peaks.first().md_Snr;
	for (int i = 1; i < lk_Peaks.size(); ++i)
		ld_MinSnr = std::min<double>(ld_MinSnr, lk_Peaks[i].md_Snr);
	
	if (ld_MinSnr < md_MinSnr)
		return lr_Result;
	
	lr_Result.md_Snr = ld_MinSnr;
	
	double ld_LightSum = 0.0;
	double ld_HeavySum = 0.0;
	
	int li_LightPeakIncludeCount = 0;
	int li_HeavyPeakIncludeCount = 0;
	int li_LightPeakExcludeCount = 0;
	int li_HeavyPeakExcludeCount = 0;
	
/*	for (int li_Isotope = 0; li_Isotope < mi_WatchIsotopesCount; ++li_Isotope)
	{
		if (ak_LightPeaksInclude.contains(li_Isotope))
		{
            if (mb_UseArea)
                ld_LightSum += ak_LightPeaksInclude[li_Isotope].md_PeakArea;
            else
                ld_LightSum += ak_LightPeaksInclude[li_Isotope].md_PeakIntensity;
			++li_LightPeakIncludeCount;
		}
		if (ak_LightPeaksExclude.contains(li_Isotope))
			++li_LightPeakExcludeCount;*/
/*		for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
		{
			if (ak_HeavyPeaksInclude.contains(k * mi_WatchIsotopesCount + li_Isotope))
			{
                if (mb_UseArea)
                    ld_HeavySum += ak_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope].md_PeakArea;
                else
                    ld_HeavySum += ak_HeavyPeaksInclude[k * mi_WatchIsotopesCount + li_Isotope].md_PeakIntensity;
				++li_HeavyPeakIncludeCount;
			}
			if (ak_HeavyPeaksExclude.contains(k * mi_WatchIsotopesCount + li_Isotope))
				++li_HeavyPeakExcludeCount;
		}*/
// 	}
	
/*	if (li_LightPeakIncludeCount == mi_WatchIsotopesCount && 
		li_HeavyPeakIncludeCount == mi_WatchIsotopesCount * mk_LabeledEnvelopeCountForPeptide[as_Peptide])
	{
		// both isotope envelopes are complete, we get a ratio!
		if (mb_CheckLightForbiddenPeaks && ak_LightPeaksExclude.contains(-1))
			// don't quantify if forbidden peak is present!
			return lr_Result;
		if (mb_CheckHeavyForbiddenPeaks)
		{
			// and don't quantify if the heavy forbidden peak is present
			for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
				if (ak_HeavyPeaksExclude.contains(-1 - k))
					return lr_Result;
		}
		lr_Result.md_AmountUnlabeled = ld_LightSum;
		lr_Result.md_AmountLabeled = ld_HeavySum;
	}
	else if (li_LightPeakIncludeCount == mi_WatchIsotopesCount && 
			 li_HeavyPeakExcludeCount == 0)
	{
		// light state only!
		if (ak_LightPeaksExclude.contains(-1))
			return lr_Result;
		lr_Result.md_AmountUnlabeled = ld_LightSum;
		lr_Result.md_AmountLabeled = 0.0;
	}
	else if (li_LightPeakExcludeCount == 0 && 
			 li_HeavyPeakIncludeCount == mi_WatchIsotopesCount * mk_LabeledEnvelopeCountForPeptide[as_Peptide])
	{
		// heavy state only
		for (int k = 0; k < mk_LabeledEnvelopeCountForPeptide[as_Peptide]; ++k)
			if (ak_HeavyPeaksExclude.contains(-1 - k))
				return lr_Result;
		lr_Result.md_AmountUnlabeled = 0.0;
		lr_Result.md_AmountLabeled = ld_HeavySum;
	}
	else
	{
		// no complete isotope envelope
		return lr_Result;
	}*/
	
	lr_Result.mk_TargetMz = ak_TargetMz;
	lr_Result.mk_ForbiddenMz = ak_ForbiddenMz;
	
	QList<double> lk_AllMz = ak_TargetMz + ak_ForbiddenMz;
	qSort(lk_AllMz);
	
	lr_Result.md_MinMz = lk_AllMz.first();
	lr_Result.md_MaxMz = lk_AllMz.last();
	
	lr_Result.mk_UnlabeledPeaks = ak_LightPeaksInclude.values();
	lr_Result.mk_LabeledPeaks = ak_HeavyPeaksInclude.values();
	lr_Result.mb_IsGood = true;
	
	return lr_Result;
}


double k_QuantifierBase::gaussian(double x, double a, double b, double c)
{
	return a * exp(-(pow(x - b, 2.0) / (2 * c * c)));
}


QHash<QString, int> k_QuantifierBase::compositionForPeptide(const QString& as_Peptide)
{
    QHash<QString, int> lk_Composition;
    lk_Composition["H"] = 2;
    lk_Composition["O"] = 1;
    for (int i = 0; i < as_Peptide.size(); ++i)
    {
        char lc_AminoAcid = as_Peptide.at(i).toAscii();
        QHash<QString, int> lk_AminoAcidComposition = mk_AminoAcidComposition[lc_AminoAcid];
        foreach (QString ls_Element, lk_AminoAcidComposition.keys())
        {
            if (!lk_Composition.contains(ls_Element))
                lk_Composition[ls_Element] = 0;
            lk_Composition[ls_Element] += lk_AminoAcidComposition[ls_Element];
        }
    }
    return lk_Composition;
}


void k_QuantifierBase::leastSquaresFit(QList<tk_DoublePair> ak_Pairs, double* ad_Factor_, double* ad_Error_)
{
    double f = 0.0;
    double e = 0.0;
    foreach (tk_DoublePair lk_Pair, ak_Pairs)
    {
        f += lk_Pair.first * lk_Pair.second;
        e += lk_Pair.first * lk_Pair.first;
    }
    f /= e;
    *ad_Factor_ = f;
    
    // determine error
    double ld_MaxTargetIntensity = 0.0;
    foreach (tk_DoublePair lk_Pair, ak_Pairs)
        ld_MaxTargetIntensity = std::max<double>(ld_MaxTargetIntensity, lk_Pair.first);
    
    double ld_Error = 0.0;

    for (int i = 0; i < ak_Pairs.size(); ++i)
    {
        tk_DoublePair lk_Pair = ak_Pairs[i];
        double ld_TargetHeight = lk_Pair.first;
        double ld_PeakHeight = lk_Pair.second;
        ld_PeakHeight /= f;
        ld_TargetHeight /= ld_MaxTargetIntensity;
        ld_PeakHeight /= ld_MaxTargetIntensity;
        ld_Error += pow(ld_PeakHeight - ld_TargetHeight, 2.0);
    }
    ld_Error /= ak_Pairs.size();
    
    *ad_Error_ = ld_Error;
}


QHash<int, int> k_QuantifierBase::matchTargetsToPeaks(QList<double> ak_PeakMz, QMultiMap<double, int> ak_Targets, double ad_MassAccuracy)
{
    // match all target m/z values simultaneously to this spectrum's peaks
    // create root bucket
    QList<double> lk_TargetMz = ak_Targets.keys();
    QList<int> lk_TargetIds = ak_Targets.values();
    
    // after parallel searching, this list will contain for 
    // each target m/z value the peak which is closest to it
    QHash<int, int> lk_PeakForTargetMz;

    int li_PeakIndex = 0;
    int li_TargetIndex = 0;
    
    while (li_TargetIndex < lk_TargetMz.size())
    {
        double ld_TargetMz = lk_TargetMz[li_TargetIndex];
        double ld_Error = ld_TargetMz / ad_MassAccuracy * 1000000.0;
        double ld_TargetMzMin = ld_TargetMz - ld_Error;
        double ld_TargetMzMax = ld_TargetMz + ld_Error;
        
        while ((li_PeakIndex < ak_PeakMz.size()) && (ak_PeakMz[li_PeakIndex] < ld_TargetMzMin))
            ++li_PeakIndex;
        
        // no more peaks, abort
        if (li_PeakIndex >= ak_PeakMz.size())
            break;
        
        double ld_PeakMz = ak_PeakMz[li_PeakIndex];
        if ((ld_PeakMz >= ld_TargetMzMin) && (ld_PeakMz <= ld_TargetMzMax))
        {
            // we found a match, but is it unambiguous?
            bool lb_Good = true;
            if (li_PeakIndex < ak_PeakMz.size() - 1)
            {
                // if there is at least one more peak, check whether it disturbs
                double ld_NextPeakMz = ak_PeakMz[li_PeakIndex + 1];
                if (ld_NextPeakMz <= ld_TargetMzMax)
                    // no, it is in the way, don't use this
                    lb_Good = false;
            }
            if (lb_Good)
                lk_PeakForTargetMz[lk_TargetIds[li_TargetIndex]] = li_PeakIndex;
            // advance target pointer
            ++li_TargetIndex;
        }
    }
    
    return lk_PeakForTargetMz;
    
    /*
    {
        // put all target m/z values into the appropriate bucket
        QList<r_Bucket> lk_NewBuckets;
        foreach (r_Bucket lr_Bucket, lk_Buckets)
        {
            if (lr_Bucket.mi_Length <= 1)
            {
                for (int i = 0; i < lr_Bucket.mk_Entries.size(); ++i)
                {
                    int li_PeakIndex = lr_Bucket.mi_Start;
                    double ld_TargetMass = lk_TargetMz[lr_Bucket.mk_Entries[i]];
                    double ld_PeakMass = ak_PeakMz[li_PeakIndex];
                    double ld_MassError = fabs(ld_PeakMass - ld_TargetMass) / ld_TargetMass * 1000000.0;
                    if (ld_MassError <= ad_MassAccuracy)
                    {
                        // the match is good, now check the adjacent peaks, if any
                        bool lb_Good = true;
                        if (li_PeakIndex > 0)
                        {
                            double ld_OtherPeakMass = ak_PeakMz[li_PeakIndex - 1];
                            double ld_OtherMassError = fabs(ld_OtherPeakMass - ld_TargetMass) / ld_TargetMass * 1000000.0;
                            if (ld_OtherMassError <= ad_MassAccuracy)
                                lb_Good = false;
                        }
                        if (li_PeakIndex < ak_PeakMz.size() - 1)
                        {
                            double ld_OtherPeakMass = ak_PeakMz[li_PeakIndex + 1];
                            double ld_OtherMassError = fabs(ld_OtherPeakMass - ld_TargetMass) / ld_TargetMass * 1000000.0;
                            if (ld_OtherMassError <= ad_MassAccuracy)
                                lb_Good = false;
                        }
                        if (lb_Good)
                            lk_PeakForTargetMz[lk_TargetIds[lr_Bucket.mk_Entries[i]]] = li_PeakIndex;
                    }
                }
            }
            else
            {
                // split this bucket and create left and right children
                int li_HalfSize = lr_Bucket.mi_Length / 2;
                r_Bucket lr_LeftChild(lr_Bucket.mi_Start, li_HalfSize);
                r_Bucket lr_RightChild(lr_Bucket.mi_Start + li_HalfSize, lr_Bucket.mi_Length - li_HalfSize);
                double ld_LeftBorder = ak_PeakMz[lr_RightChild.mi_Start - 1];
                double ld_RightBorder = ak_PeakMz[lr_RightChild.mi_Start];
                double ld_Razor = (ld_LeftBorder + ld_RightBorder) * 0.5;
                
                // now determine which target m/z entries go into the 
                // left child and which go into the right child
                
                // TODO: this could be sped up because everything is sorted
                for (int i = 0; i < lr_Bucket.mk_Entries.size(); ++i)
                {
                    if (lk_TargetMz[lr_Bucket.mk_Entries[i]] > ld_Razor)
                        lr_RightChild.mk_Entries.push_back(lr_Bucket.mk_Entries[i]);
                    else
                        lr_LeftChild.mk_Entries.push_back(lr_Bucket.mk_Entries[i]);
                }
                
                if (!lr_LeftChild.mk_Entries.empty())
                    lk_NewBuckets.push_back(lr_LeftChild);
                if (!lr_RightChild.mk_Entries.empty())
                    lk_NewBuckets.push_back(lr_RightChild);
            }
        }
        lk_Buckets = lk_NewBuckets;
    }    
    return lk_PeakForTargetMz;
    */
}


QHash<int, int> k_QuantifierBase::extractMatches(QSet<int> ak_Ids, QList<int> ak_TargetIdsSorted, QHash<int, int> ak_Matches)
{
    QHash<int, int> lk_Result;
    
    // ak_Matches contains a peak id for a target id
    foreach (int li_Id, ak_Matches.keys())
    {
        int li_TranslatedId = ak_TargetIdsSorted[li_Id];
        if (ak_Ids.contains(li_TranslatedId))
            lk_Result[li_TranslatedId] = ak_Matches[li_Id];
    }
    
    return lk_Result;
}


void k_QuantifierBase::parseLabel()
{
    QStringList lk_Tokens = tokenize(ms_Label.toUpper());

    int li_State = 0;
    QSet<QString> lk_AminoAcidScope;
    
    // this hash contains special environmental conditions for every amino acid
    // if an amino acid is not contained, it comes from natural conditions!
    QHash<QString, tk_ArtificialEnvironment> lk_EnvironmentForAminoAcid;

    while (!lk_Tokens.empty())
    {
        QVariant::Type le_TokenType;
        QString ls_Token = fetchNextToken(&lk_Tokens, &le_TokenType);
        if (li_State == 0)
        {
            if (le_TokenType == QVariant::String)
            {
                // we found a scope-defining amino acid prefix
                lk_AminoAcidScope << ls_Token;
            }
            else
                li_State = 1;
        }
        if (li_State == 1)
        {
            if (le_TokenType != QVariant::Int)
            {
                printf("Error: Isotope number expected, but got %s.\n", ls_Token.toStdString().c_str());
                exit(1);
            }
            // fetch element
            if (lk_Tokens.empty())
            {
                printf("Error: Element expected, but string is at end.\n");
                exit(1);
            }
            
            int li_Isotope = ls_Token.toInt();
            
            ls_Token = fetchNextToken(&lk_Tokens, &le_TokenType);
            if (le_TokenType != QVariant::String)
            {
                printf("Error: Element expected, but got %s.\n", ls_Token.toStdString().c_str());
                exit(1);
            }
            
            QString ls_Element = ls_Token;
            double ld_Efficiency = 1.0;
            
            if (peekNextToken(lk_Tokens) == QVariant::Double)
            {
                // parse labeling efficiency
                ls_Token = fetchNextToken(&lk_Tokens, &le_TokenType);
                ld_Efficiency = ls_Token.toFloat();
                if (!(ls_Element == "C" || ls_Element == "N"))
                {
                    printf("Error: A labeling efficiency can only be specified for C and N.\n");
                    exit(1);
                }
            }
            
            if (lk_AminoAcidScope.empty() || lk_AminoAcidScope.contains("X"))
            {
                lk_AminoAcidScope = QSet<QString>();
                QString ls_AllAminoAcids = "GASPVTCLINDQKEMHFRYW";
                for (int i = 0; i < ls_AllAminoAcids.length(); ++i)
                    lk_AminoAcidScope << ls_AllAminoAcids.mid(i, 1);
            }
            foreach (QString ls_AminoAcid, lk_AminoAcidScope)
            {
                if (!lk_EnvironmentForAminoAcid.contains(ls_AminoAcid))
                    lk_EnvironmentForAminoAcid[ls_AminoAcid] = tk_ArtificialEnvironment();
                lk_EnvironmentForAminoAcid[ls_AminoAcid][ls_Element] = r_IsotopeAbundance(li_Isotope - mk_IsotopeEnvelope.mk_BaseIsotope[ls_Element], ld_Efficiency);
            }
            
            if (peekNextToken(lk_Tokens) != QVariant::Int)
            {
                lk_AminoAcidScope.clear();
                li_State = 0;
            }
        }
    }
    
    QMultiMap<QString, QString> lk_AminoAcidForDescription;
    
    foreach (QString ls_AminoAcid, lk_EnvironmentForAminoAcid.keys())
    {
//         printf("AE for %s:\n", ls_AminoAcid.toStdString().c_str());
        tk_ArtificialEnvironment lk_Environment = lk_EnvironmentForAminoAcid[ls_AminoAcid];
        QHash<QString, QList<double> > lk_Abundances;
        foreach (QString ls_Element, lk_Environment.keys())
        {
            r_IsotopeAbundance lr_IsotopeAbundance = lk_Environment[ls_Element];
//             printf("%d%s at %1.2f\n", lr_IsotopeAbundance.mi_Isotope, ls_Element.toStdString().c_str(), lr_IsotopeAbundance.mf_Efficiency);
            QList<double> lk_AbundanceNumbers;
            lk_AbundanceNumbers << 0.0 << 0.0 << 0.0 << 0.0 << 0.0;
            lk_AbundanceNumbers[lr_IsotopeAbundance.mi_Isotope] = lr_IsotopeAbundance.mf_Efficiency;
            if (lr_IsotopeAbundance.mf_Efficiency < 1.0)
                lk_AbundanceNumbers[1 - lr_IsotopeAbundance.mi_Isotope] = 1.0 - lr_IsotopeAbundance.mf_Efficiency;
            lk_Abundances[ls_Element] = lk_AbundanceNumbers;
        }
        mk_HeavyIsotopeEnvelopeForAminoAcid[ls_AminoAcid] = k_IsotopeEnvelope(lk_Abundances);
        
        // make description
        QStringList lk_ElementList = lk_Abundances.keys();
        qSort(lk_ElementList.begin(), lk_ElementList.end());
        QStringList lk_Description;
        foreach (QString ls_Element, lk_ElementList)
        {
            for (int i = 0; i < lk_Abundances[ls_Element].size(); ++i)
            {
                double ld_Abundance = lk_Abundances[ls_Element][i];
                if (ld_Abundance > 0.0)
                    lk_Description << QString("%1% %2%3")
                        .arg(ld_Abundance * 100, 0, 'f', 1)
                        .arg(i + mk_IsotopeEnvelope.mk_BaseIsotope[ls_Element])
                        .arg(ls_Element);
            }
        }
        lk_AminoAcidForDescription.insert(lk_Description.join(", "), ls_AminoAcid);
    }
/*    if (!mb_Quiet)
    {
        printf("Label composition:\n");
        foreach (QString ls_Description, lk_AminoAcidForDescription.uniqueKeys())
        {
            QStringList lk_AminoAcids = lk_AminoAcidForDescription.values(ls_Description);
            if (lk_AminoAcids.size() == 20)
                lk_AminoAcids = QStringList() << "all amino acids";
            printf("%s: %s\n", lk_AminoAcids.join(", ").toStdString().c_str(), ls_Description.toStdString().c_str());
        }
    }*/
}


QStringList k_QuantifierBase::tokenize(QString as_String)
{
    QStringList lk_Result;
    bool lb_Append = true;
    for (int i = 0; i < as_String.length(); ++i)
    {
        QChar lk_Char = as_String.at(i);
        if (lk_Char.isSpace())
            lb_Append = false;
        if (lk_Char.isDigit() || lk_Char == QChar('.'))
        {
            if (lk_Result.empty())
            {
                lk_Result << QString(lk_Char);
                lb_Append = true;
            }
            else
            {
                if (lb_Append & (lk_Result.last().at(0).isDigit() || lk_Result.last().at(0) == QChar('.')))
                    lk_Result.last() += lk_Char;
                else
                {
                    lk_Result << QString(lk_Char);
                    lb_Append = true;
                }
            }
        }
        else if (lk_Char.isLetter())
        {
            lk_Result << QString(lk_Char);
            lb_Append = true;
        }
    }
    return lk_Result;
}


QString k_QuantifierBase::fetchNextToken(QStringList* ak_StringList_, QVariant::Type* ae_Type_)
{
    QVariant::Type le_TokenType = QVariant::String;
    QString ls_Token = ak_StringList_->takeFirst();
    bool lb_Ok = false;
    ls_Token.toFloat(&lb_Ok);
    if (lb_Ok)
    {
        if (ls_Token.contains("."))
            le_TokenType = QVariant::Double;
        else
            le_TokenType = QVariant::Int;
    }
    *ae_Type_ = le_TokenType;
    return ls_Token;
}


QVariant::Type k_QuantifierBase::peekNextToken(QStringList ak_StringList)
{
    if (ak_StringList.empty())
        return QVariant::Invalid;
    QVariant::Type le_TokenType = QVariant::String;
    QString ls_Token = ak_StringList.first();
    bool lb_Ok = false;
    ls_Token.toFloat(&lb_Ok);
    if (lb_Ok)
    {
        if (ls_Token.contains("."))
            le_TokenType = QVariant::Double;
        else
            le_TokenType = QVariant::Int;
    }
    return le_TokenType;
}


tk_IsotopeEnvelope k_QuantifierBase::lightEnvelopeForPeptide(QString as_Peptide)
{
    QHash<QString, int> lk_Composition = compositionForPeptide(as_Peptide);
    return mk_IsotopeEnvelope.isotopeEnvelopeForComposition(lk_Composition);
}


tk_IsotopeEnvelope k_QuantifierBase::heavyEnvelopeForPeptide(QString as_Peptide)
{
    QHash<QString, int> lk_Composition = compositionForPeptide(as_Peptide);
    
    QHash<QString, QHash<QString, int> > lk_CompositionForAminoAcid;
    for (int i = 0; i < as_Peptide.length(); ++i)
    {
        QString ls_AminoAcid = as_Peptide.mid(i, 1);
        QString ls_Id = QString();
        if (mk_HeavyIsotopeEnvelopeForAminoAcid.contains(ls_AminoAcid))
            ls_Id = ls_AminoAcid;
        QHash<QString, int> lk_ThisComposition = compositionForPeptide(ls_AminoAcid);
        foreach (QString ls_Element, lk_ThisComposition.keys())
        {
            if (!lk_CompositionForAminoAcid[ls_Id].contains(ls_Element))
                lk_CompositionForAminoAcid[ls_Id][ls_Element] = 0;
            lk_CompositionForAminoAcid[ls_Id][ls_Element] += lk_ThisComposition[ls_Element];
        }
    }
    
    tk_IsotopeEnvelope lr_HeavyMass;
    bool lb_First = true;
    foreach (QString ls_Id, lk_CompositionForAminoAcid.keys())
    {
        k_IsotopeEnvelope* lk_UseEnvelope_ = &mk_IsotopeEnvelope;
        if (!ls_Id.isEmpty())
            lk_UseEnvelope_ = &mk_HeavyIsotopeEnvelopeForAminoAcid[ls_Id];
        if (lb_First)
            lr_HeavyMass = lk_UseEnvelope_->isotopeEnvelopeForComposition(lk_CompositionForAminoAcid[ls_Id]);
        else
            lr_HeavyMass = mk_IsotopeEnvelope.add(lr_HeavyMass, lk_UseEnvelope_->isotopeEnvelopeForComposition(lk_CompositionForAminoAcid[ls_Id]));
        lb_First = false;
    }
    
    return lr_HeavyMass;
}
