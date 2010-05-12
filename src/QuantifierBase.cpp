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
    , md_AbsenceMassAccuracyFactor(DEFAULT_ABSENCE_MASS_ACCURACY_FACTOR)
    , md_RequireAbundance(DEFAULT_REQUIRE_ABUNDANCE)
    , md_ConsiderAbundance(DEFAULT_CONSIDER_ABUNDANCE)
    , md_MaxFitError(DEFAULT_MAX_FIT_ERROR)
    , mi_FixedIsotopePeakCount(DEFAULT_FIXED_ISOTOPE_PEAK_COUNT)
    , mb_CheckForbiddenPeak(DEFAULT_CHECK_FORBIDDEN_PEAK)
    , mb_CheckOverlappingPeaks(DEFAULT_CHECK_OVERLAPPING_PEAKS)
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
    if (md_AbsenceMassAccuracyFactor < 0.0)
    {
        printf("Error: absence mass accuracy factor must not be less than zero.\n");
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
    if (md_MaxFitError < 0.0)
    {
        printf("Error: max fit error must not be less than zero.\n");
        exit(1);
    }
    if (mi_FixedIsotopePeakCount < 1)
    {
        printf("Error: Fixed isotope peak count must not be less than 1.\n");
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
    
    foreach (char lc_AminoAcid, mk_AminoAcidComposition.keys())
        mk_LightIsotopeEnvelopeForAminoAcid[QString("%1").arg(lc_AminoAcid)] = mk_IsotopeEnvelope.isotopeEnvelopeForComposition(compositionForPeptide(QString("%1").arg(lc_AminoAcid)));
    
    if (mk_pCsvDevice)
        mk_pCsvStream = QSharedPointer<QTextStream>(new QTextStream(mk_pCsvDevice.data()));
    if (mk_pXhtmlDevice)
        mk_pXhtmlStream = QSharedPointer<QTextStream>(new QTextStream(mk_pXhtmlDevice.data()));
    
    parseLabel();
}


k_QuantifierBase::~k_QuantifierBase()
{
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
    if (ak_Arguments.empty())
        printUsageAndExit();
    
    if (!ak_Arguments.empty() && ak_Arguments.first() == "--version")
    {
        printf("%s %s\n", ms_ProgramName.toStdString().c_str(), gs_Version.toStdString().c_str());
        exit(0);
    }
    
    if (!ak_Arguments.empty() && ak_Arguments.first() == "--help")
        printUsageAndExit();
    
    QFile* lk_StdOut_ = new QFile();
    lk_StdOut_->open(stdout, QIODevice::WriteOnly);
    mk_pCsvDevice = QSharedPointer<QIODevice>(lk_StdOut_);

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
                mk_pCsvDevice = QSharedPointer<QIODevice>();
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
            mk_pCsvDevice = QSharedPointer<QIODevice>(lk_File_);
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
            mk_pXhtmlDevice = QSharedPointer<QIODevice>(lk_File_);
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
        else if (mk_Parameters.contains(r_Parameter::MassAccuracy) && ls_Key == "--absenceMassAccuracyFactor")
            md_AbsenceMassAccuracyFactor = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::RequireAbundance) && ls_Key == "--requireAbundance")
            md_RequireAbundance = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::ConsiderAbundance) && ls_Key == "--considerAbundance")
            md_ConsiderAbundance = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MaxFitError) && ls_Key == "--maxFitError")
            md_MaxFitError = stringToDouble(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::MaxFitError) && ls_Key == "--isotopePeaks")
            mi_FixedIsotopePeakCount = stringToInt(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::CheckForbiddenPeak) && ls_Key == "--checkForbiddenPeak")
            mb_CheckForbiddenPeak = stringToBool(ak_Arguments.takeFirst());
        else if (mk_Parameters.contains(r_Parameter::CheckOverlappingPeaks) && ls_Key == "--checkOverlappingPeaks")
            mb_CheckOverlappingPeaks = stringToBool(ak_Arguments.takeFirst());
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
}


double k_QuantifierBase::calculatePeptideMass(QString as_Peptide, int ai_Charge)
{
    double ld_Mass = md_WaterMass;
    for (int i = 0; i < as_Peptide.length(); ++i)
        ld_Mass += mk_AminoAcidWeight[as_Peptide.at(i).toAscii()];
    return (ld_Mass + md_HydrogenMass * ai_Charge) / ai_Charge;
}


QString k_QuantifierBase::renderScanAsSvg(r_Scan& ar_Scan, r_ScanQuantitationResult ar_QuantitationResult, QList<r_Peak>& ak_AllPeaks, QHash<int, int>& ak_Matches)
{
    double ld_Ratio = 4.0;
    double ld_Width = 950.0;
    double ld_Height = ld_Width / ld_Ratio;
    double ld_BorderTop = 4.0;
    double ld_BorderRight = 16.0;
    double ld_BorderBottom = 20.0;
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
    double ymax = 0.0;
/*    foreach (QString ls_Title, ar_QuantitationResult.mk_UnlabeledProfileScale.keys())
    {
        double ld_Scale = ar_QuantitationResult.mk_UnlabeledProfileScale[ls_Title];
        if (ld_Scale > 0.0)
        {
            tk_IsotopeEnvelope lk_IsotopeEnvelope = mk_RenderIsotopeEnvelopeForPeptideWeight[ar_QuantitationResult.ms_Peptide + "-0"][ls_Title];
            for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
            {
                double ld_Abundance = lk_IsotopeEnvelope[i].first * ld_Scale;
                ymax = std::max<double>(ymax, ld_Abundance);
            }
        }
    }
    foreach (QString ls_Title, ar_QuantitationResult.mk_LabeledProfileScale.keys())
    {
        double ld_Scale = ar_QuantitationResult.mk_LabeledProfileScale[ls_Title];
        if (ld_Scale > 0.0)
        {
            tk_IsotopeEnvelope lk_IsotopeEnvelope = mk_RenderIsotopeEnvelopeForPeptideWeight[ar_QuantitationResult.ms_Peptide + "-1"][ls_Title];
            for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
            {
                double ld_Abundance = lk_IsotopeEnvelope[i].first * ld_Scale;
                ymax = std::max<double>(ymax, ld_Abundance);
            }
        }
    }*/
    for (int li_Weight = 0; li_Weight < 2; ++li_Weight)
    {
        QString ls_Key = QString("%1-%2-%3").arg(ar_QuantitationResult.ms_Peptide).arg(ar_QuantitationResult.mi_Charge).arg(li_Weight);
        foreach (r_EnvelopePeaks lr_Peaks, mk_TargetsForPeptideChargeWeight[ls_Key])
        {
            foreach (int li_Id, lr_Peaks.mk_RequiredIds | lr_Peaks.mk_ConsideredIds)
            {
                if (ak_Matches.contains(li_Id))
                {
                    r_Peak& lr_Peak = ak_AllPeaks[ak_Matches[li_Id]];
                    ymax = std::max<double>(ymax, lr_Peak.md_PeakIntensity);
                }
            }
        }
    }
    
    // add a little space above
    ymax *= 1.1;
            
    // ymax remains unscaled!
    double ymaxScaled = this->scale(1.0);
    
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

        tk_IsotopeEnvelope lk_IsotopeEnvelope;
        bool lb_DrawnOnePoint;
        double ld_LastMz;
        
        double ld_PeptideBaseMass = mk_RenderBaseMassForPeptide[ls_Peptide];
        
        if (mb_UseIsotopeEnvelopes)
        {
            foreach (QString ls_WeightTitle, ar_QuantitationResult.mk_ProfileScale.keys())
            {
                double ld_ProfileScale = ar_QuantitationResult.mk_ProfileScale[ls_WeightTitle];
                double lb_Good = ar_QuantitationResult.mk_Good[ls_WeightTitle];
                if (ld_ProfileScale > 0.0)
                {
                    lk_IsotopeEnvelope = mk_RenderIsotopeEnvelopeForPeptideWeightTitle[ls_Peptide + "-" + ls_WeightTitle];
                    
                    lb_DrawnOnePoint = false;
                    ld_LastMz = 0.0;
                    QPainterPath lk_Path;
                    for (int i = 0; i < lk_IsotopeEnvelope.size(); ++i)
                    {
                        double ld_Abundance = lk_IsotopeEnvelope[i].first * ld_ProfileScale;
                        double ld_Mz = (ld_PeptideBaseMass + lk_IsotopeEnvelope[i].second + md_HydrogenMass * li_Charge) / li_Charge;
                        double x = (ld_Mz - xmin) * dx + x0;
                        double y = this->scale(ld_Abundance / ymax) / ymaxScaled;
                        if (y > 0.001)
                        {
                            y = y * dy + y0;
                            if (!lb_DrawnOnePoint)
                            {
                                lk_Path.moveTo(QPointF(((ld_Mz - 0.5 / li_Charge) - xmin) * dx + x0, y0));
                                lb_DrawnOnePoint = true;
                            }
                            lk_Path.lineTo(QPointF(x, y));
                            ld_LastMz = ld_Mz;
                        }
                    }
                    lk_Path.lineTo(QPointF(((ld_LastMz + 0.5 / li_Charge) - xmin) * dx + x0, y0));
                    lk_Pen.setWidthF(1.0);
                    lk_Pen.setJoinStyle(Qt::BevelJoin);
                    lk_Pen.setColor(QColor(TANGO_ALUMINIUM_3));
                    if (lb_Good)
                    {
                        lk_Pen.setStyle(Qt::SolidLine);
                        lk_Painter.setBrush(QBrush(QColor(TANGO_ALUMINIUM_0)));
                    }
                    else
                    {
                        lk_Pen.setStyle(Qt::DashLine);
                        lk_Painter.setBrush(Qt::NoBrush);
                    }
                    lk_Painter.setPen(lk_Pen);
                    
                    lk_Painter.drawPath(lk_Path);
                }
            }
        }

        lk_Painter.setBrush(Qt::NoBrush);
        // draw spectrum
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::RoundJoin);
        lk_Pen.setColor(QColor(TANGO_CHAMELEON_2));
        lk_Pen.setStyle(Qt::DashLine);
        lk_Painter.setPen(lk_Pen);
        
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
        
        // draw x ticks
        
/*        // xt0 and xts are x tick start and step
        double xt0 = xmin;
        double xts = 1.0;
        // determine an appropriate step value
        
        double xt = xt0;
        while (xt < xmax)
        {
            lk_Painter.drawLine(QPointF((xt - xmin) * dx + x0, y0), QPointF((xt - xmin) * dx + x0, y0 + 4.0));
            ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana; text-anchor: middle;'>%3</text>\n").arg((xt - xmin) * dx + x0).arg(y0 + 14.0).arg((int)round(xt));
            xt += xts;
        }*/
        
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::RoundJoin);
        lk_Pen.setColor(QColor(192, 192, 192));
        lk_Pen.setStyle(Qt::SolidLine);
        lk_Painter.setPen(lk_Pen);
        
        QVector<QPointF> lk_Points;
        for (int i = li_Start; i <= li_End; ++i)
            lk_Points.append(QPointF((ar_Scan.mr_Spectrum.md_MzValues_[i] - xmin) * dx + x0, 
                                    this->scale(ar_Scan.mr_Spectrum.md_IntensityValues_[i] / ymax) / ymaxScaled * dy + y0));
        lk_Painter.drawPolyline(QPolygonF(lk_Points));
        
        lk_Pen.setWidthF(1.0);
        lk_Pen.setJoinStyle(Qt::RoundJoin);
        lk_Pen.setColor(QColor(TANGO_CHAMELEON_2));
        lk_Pen.setStyle(Qt::SolidLine);
        lk_Painter.setPen(lk_Pen);
        
        for (int li_Weight = 0; li_Weight < 2; ++li_Weight)
        {
            QString ls_Key = QString("%1-%2-%3").arg(ar_QuantitationResult.ms_Peptide).arg(ar_QuantitationResult.mi_Charge).arg(li_Weight);
            foreach (r_EnvelopePeaks lr_Peaks, mk_TargetsForPeptideChargeWeight[ls_Key])
            {
                double x = (lr_Peaks.md_BaseMz - xmin) * dx + x0;
                //lk_Painter.drawLine(x, y0), QPointF(x, y0 + 4.0));
                ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana;'>%3</text>\n").arg(x - 4.0).arg(y0 + 14.0).arg(QString("%1 (%2)").arg(lr_Peaks.md_BaseMz, 1, 'f', 3).arg(lr_Peaks.ms_Title));
                foreach (int li_Id, lr_Peaks.mk_RequiredIds | lr_Peaks.mk_ConsideredIds)
                {
                    if (ak_Matches.contains(li_Id))
                    {
                        if (lr_Peaks.mk_RequiredIds.contains(li_Id))
                            lk_Pen.setStyle(Qt::SolidLine);
                        else
                            lk_Pen.setStyle(Qt::DotLine);
                        lk_Pen.setColor(QColor(TANGO_SKY_BLUE_1));
                        lk_Painter.setPen(lk_Pen);
                        
                        r_Peak& lr_Peak = ak_AllPeaks[ak_Matches[li_Id]];
                        double x = (lr_Peak.md_PeakMz - xmin) * dx + x0;
                        double y = this->scale(lr_Peak.md_PeakIntensity / ymax) / ymaxScaled * dy + y0;
                        lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y));
                        lk_Painter.drawArc(QRectF(QPointF(x, y) - QPointF(3.0, 3.0), QSize(6.0, 6.0)), 0, 5760);
                    }
                    else
                    {
                        if (lr_Peaks.mk_RequiredIds.contains(li_Id))
                        {
                            lk_Pen.setStyle(Qt::SolidLine);
                            lk_Pen.setColor(QColor(TANGO_SCARLET_RED_2));
                            lk_Painter.setPen(lk_Pen);
                            
                            double x = (mk_TargetMzAndIntensity[li_Id].first - xmin) * dx + x0;
                            lk_Painter.drawLine(QPointF(x, y0), QPointF(x, y0 + dy));
                        }
                    }
                }
            }
        }
        
    }
    lk_Buffer.close();
    lk_Buffer.open(QBuffer::ReadOnly);
    lk_Buffer.seek(0);
    QString ls_Result = QString(lk_Buffer.readAll());
    lk_Buffer.close();
/*    foreach (r_Peak lr_Peak, (ar_QuantitationResult.mk_UnlabeledPeaks + ar_QuantitationResult.mk_LabeledPeaks))
    {
        double ld_X, ld_Y;
        ld_X = (lr_Peak.md_PeakMz - xmin) * dx + x0;
        ld_Y = this->scale(lr_Peak.md_PeakIntensity / ymax) / ymaxScaled * dy + y0;
        char lc_Number_[1024];
        lc_Number_[0] = 0;
        // :TODO: now this always uses peak height, make this an option (area/intensity)
        sprintf(lc_Number_, "%.2g", lr_Peak.md_PeakIntensity);
        if (strlen(lc_Number_) > 0)
            ls_Labels += QString("<text x='%1' y='%2' style='font-size:10px; font-family:Verdana; fill: #aaa;' transform='rotate(90 %1 %2)'>%3</text>\n").arg(ld_X + 3.0).arg(ld_Y + 3.0).arg(lc_Number_);
    }*/
    ls_Result.replace("</svg>", ls_Labels + "</svg>");
    return ls_Result;
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
    
    // MAIN SECTION
    
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
    if (mk_Parameters.contains(r_Parameter::UseIsotopeEnvelopes))
    {
        printf("  --useIsotopeEnvelopes [yes|no] (default: %s)\n", DEFAULT_USE_ISOTOPE_ENVELOPES ? "yes" : "no");
        printf("      Specify whether amounts should be estimated via isotope envelope fitting\n");
        printf("      or whether a constant number of peaks should be used instead.\n");
    }

    // PEAK PICKING SECTION

    if (mk_Parameters.contains(r_Parameter::ScanType))
    {
        printf("  --scanType [full|sim|all] (default: all)\n");
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
    if (mk_Parameters.contains(r_Parameter::CheckForbiddenPeak))
    {
        printf("  --checkForbiddenPeak [yes|no] (default: %s)\n", DEFAULT_CHECK_FORBIDDEN_PEAK ? "yes" : "no");
        printf("      Specify whether the forbidden peak is required to be absent.\n");
    }
    if (mk_Parameters.contains(r_Parameter::CheckOverlappingPeaks))
    {
        printf("  --checkOverlappingPeaks [yes|no] (default: %s)\n", DEFAULT_CHECK_OVERLAPPING_PEAKS ? "yes" : "no");
        printf("      Specify whether overlapping peaks should be checked for.\n");
    }
    if (mk_Parameters.contains(r_Parameter::AbsenceMassAccuracyFactor))
    {
        printf("  --absenceMassAccuracyFactor <float> (default: %1.1f)\n", DEFAULT_ABSENCE_MASS_ACCURACY_FACTOR);
        printf("      For the forbidden peak, a higher mass accuracy may be specified\n");
        printf("      this factor. The specified mass accuracy is multiplied with this factor\n");
        printf("      when checking for peak absence.\n");
    }

    // ISOTOPE ENVELOPE FITTING SECTION

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
    
    // FIXED ISOTOPE PEAK COUNT SECTION

    if (mk_Parameters.contains(r_Parameter::FixedIsotopePeakCount))
    {
        printf("  --isotopePeaks <int> (default: %d)\n", DEFAULT_FIXED_ISOTOPE_PEAK_COUNT);
        printf("      Specify how many isotope peaks should be added for amount estimation,\n");
        printf("      if --useIsotopeEnvelopes is set to 'no'.\n");
    }
    
    // OUTPUT FILES SECTION
    
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
    
    // TWEAKS SECTION
    
    if (mk_Parameters.contains(r_Parameter::LogScale))
    {
        printf("  --logScale [yes|no] (default: %s)\n", DEFAULT_LOG_SCALE ? "yes" : "no");
        printf("      Use logarithmic scale in XHTML spectra.\n");
    }
    
    // MISC. FLAGS SECTION
    
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


void k_QuantifierBase::leastSquaresFit(QList<tk_DoublePair> ak_Pairs, double* ad_Factor_, 
                                       QList<double>* ak_Errors_)
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
    
    for (int i = 0; i < ak_Pairs.size(); ++i)
    {
        tk_DoublePair lk_Pair = ak_Pairs[i];
        double ld_TargetHeight = lk_Pair.first;
        double ld_PeakHeight = lk_Pair.second;
        
        // multiply peaks with fit factor
        ld_PeakHeight /= f;

        // normalize so that max target intensity is 1.0
        ld_TargetHeight /= ld_MaxTargetIntensity;
        ld_PeakHeight /= ld_MaxTargetIntensity;
        
        double ld_ThisError = fabs(ld_PeakHeight - ld_TargetHeight);
        (*ak_Errors_) << ld_ThisError;
    }
}


QHash<int, int> k_QuantifierBase::matchTargetsToPeaks(QList<double> ak_PeakMz, QMultiMap<double, int> ak_TargetMzMin, QHash<int, tk_DoublePair> ak_TargetMzAndIntensity, QSet<int> ak_ForbiddenIds)
{
    // match all target m/z values simultaneously to this spectrum's peaks
    // create root bucket
    QList<double> lk_TargetMzMin = ak_TargetMzMin.keys();
    QList<int> lk_TargetIds = ak_TargetMzMin.values();
    
    // after the search, this list will contain for 
    // each target m/z value the peak which is closest to it
    QHash<int, int> lk_PeakForTargetMz;
    
    if (lk_TargetMzMin.empty() || ak_PeakMz.empty())
        return lk_PeakForTargetMz;
    
    int li_PeakIndex = 0;
    int li_TargetIndex = 0;
    double ld_TargetMzMin = lk_TargetMzMin[li_TargetIndex];
    double ld_TargetMzMax = 2.0 * ak_TargetMzAndIntensity[lk_TargetIds[li_TargetIndex]].first - ld_TargetMzMin;
    
    while (true)
    {
        // advance both peak pointers while we're not within the target range
        while (ak_PeakMz[li_PeakIndex] < ld_TargetMzMin)
        {
            ++li_PeakIndex;
            if (li_PeakIndex >= ak_PeakMz.size())
                break;
        }
        if (li_PeakIndex >= ak_PeakMz.size())
            break;
        
        if (ak_PeakMz[li_PeakIndex] <= ld_TargetMzMax)
        {
            // we're within target range now
            // make sure that there is no other peak within the target range
            if ((ak_ForbiddenIds.contains(li_TargetIndex)) || ((li_PeakIndex >= ak_PeakMz.size() - 1) || (ak_PeakMz[li_PeakIndex + 1] > ld_TargetMzMax)))
            {
                // assign the peak, because
                // 1. it's a forbidden target, and it doesn't matter how many peaks match, or
                // 2. it's the last peak, or
                // 3. there's only one peak within the target range
                lk_PeakForTargetMz[li_TargetIndex] = li_PeakIndex;
            }
        }
        
        // advance target pointer
        ++li_TargetIndex;
        if (li_TargetIndex >= lk_TargetMzMin.size())
            break;
        
        ld_TargetMzMin = lk_TargetMzMin[li_TargetIndex];
        ld_TargetMzMax = 2.0 * ak_TargetMzAndIntensity[lk_TargetIds[li_TargetIndex]].first - ld_TargetMzMin;
    }
    
    return lk_PeakForTargetMz;
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
    bool lb_AminoAcidScopeNegated = false;
    QString ls_LastAminoAcid;
    
    QString ls_AllAminoAcids = "GASPVTCLINDQKEMHFRYW";
    QSet<QString> lk_AllAminoAcids;
    for (int i = 0; i < ls_AllAminoAcids.length(); ++i)
        lk_AllAminoAcids << ls_AllAminoAcids.mid(i, 1);
    
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
                if (ls_Token == "*")
                {
                    if (ls_LastAminoAcid.isEmpty())
                    {
                        printf("Error: A star (*) must follow immediately after an amino acid.\n");
                        exit(1);
                    }
                    mk_StarAminoAcids << ls_LastAminoAcid;
                }
                else if (ls_Token == "^")
                {
                    if (!lk_AminoAcidScope.empty())
                    {
                        printf("Error: An amino acid scope negator (^) is only allowed as the first symbol.\n");
                        exit(1);
                    }
                    lb_AminoAcidScopeNegated = true;
                }
                else
                {
                    ls_LastAminoAcid = ls_Token;
                    lk_AminoAcidScope << ls_Token;
                }
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
                for (int i = 0; i < ls_AllAminoAcids.length(); ++i)
                    lk_AminoAcidScope << ls_AllAminoAcids.mid(i, 1);
            }
            QSet<QString> lk_UseAminoAcidScope = lk_AminoAcidScope;
            if (lb_AminoAcidScopeNegated)
                lk_UseAminoAcidScope = lk_AllAminoAcids - lk_AminoAcidScope;
            foreach (QString ls_AminoAcid, lk_UseAminoAcidScope)
            {
                if (!lk_EnvironmentForAminoAcid.contains(ls_AminoAcid))
                    lk_EnvironmentForAminoAcid[ls_AminoAcid] = tk_ArtificialEnvironment();
                lk_EnvironmentForAminoAcid[ls_AminoAcid][ls_Element] = r_IsotopeAbundance(li_Isotope - mk_IsotopeEnvelope.mk_BaseIsotope[ls_Element], ld_Efficiency);
//                 mk_NominalMassShiftForAminoAcid[ls_AminoAcid]
            }
            
            if (peekNextToken(lk_Tokens) != QVariant::Int)
            {
                // reset everything
                lk_AminoAcidScope.clear();
                li_State = 0;
                lb_AminoAcidScopeNegated = false;
            }
        }
    }
    
    foreach (QString ls_AminoAcid, lk_EnvironmentForAminoAcid.keys())
    {
//         printf("AE for %s:\n", ls_AminoAcid.toStdString().c_str());
        tk_ArtificialEnvironment lk_Environment = lk_EnvironmentForAminoAcid[ls_AminoAcid];
        QHash<QString, QList<double> > lk_Abundances;
        QHash<QString, QList<double> > lk_AbundancesFull;
        double ld_MassShift = 0.0;
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
            int li_AtomCount = 0;
            if (mk_AminoAcidComposition[ls_AminoAcid.at(0).toAscii()].contains(ls_Element))
                li_AtomCount = mk_AminoAcidComposition[ls_AminoAcid.at(0).toAscii()][ls_Element];
            if (li_AtomCount > 0)
                ld_MassShift += mk_IsotopeEnvelope.mk_ElementIsotopeMassShift[ls_Element][lr_IsotopeAbundance.mi_Isotope] * li_AtomCount;
        }
        mk_HeavyIsotopeEnvelopeForAminoAcid[ls_AminoAcid] = k_IsotopeEnvelope(lk_Abundances).isotopeEnvelopeForComposition(compositionForPeptide(ls_AminoAcid));
        mk_NominalMassShiftForAminoAcid[ls_AminoAcid] = ld_MassShift;
        
        // make description
        QStringList lk_ElementList = lk_Abundances.keys();
        qSort(lk_ElementList.begin(), lk_ElementList.end());
        QStringList lk_Description;
        foreach (QString ls_Element, lk_ElementList)
        {
            for (int i = 0; i < lk_Abundances[ls_Element].size(); ++i)
            {
                double ld_Abundance = lk_Abundances[ls_Element][i];
                if (ld_Abundance > 0.0 && i > 0)
                    lk_Description << QString("%1%2%3")
                        .arg(i + mk_IsotopeEnvelope.mk_BaseIsotope[ls_Element])
                        .arg(ls_Element)
                        .arg(ld_Abundance != 1.0 ? QString(" (%1%)").arg(ld_Abundance * 100, 0, 'f', 1) : "");
            }
        }
        mk_AminoAcidForDescription.insert(lk_Description.join(", "), ls_AminoAcid);
    }
    if (!mb_Quiet)
    {
        printf("Label composition:\n");
        foreach (QString ls_Description, mk_AminoAcidForDescription.uniqueKeys())
        {
            QStringList lk_AminoAcids = mk_AminoAcidForDescription.values(ls_Description);
            qSort(lk_AminoAcids);
            if (lk_AminoAcids.size() == 20)
                lk_AminoAcids = QStringList() << "all amino acids";
            QString ls_Scope = lk_AminoAcids.join("");
            if (ls_Scope.size() > 1)
                ls_Scope = "[" + ls_Scope + "]";
            printf("%s: %s\n", ls_Scope.toStdString().c_str(), ls_Description.toStdString().c_str());
        }
    }
    if ((!mk_StarAminoAcids.empty()) && mb_UseIsotopeEnvelopes)
    {
        printf("Warning: You are using a label containing the star symbol (*) in conjunction\n");
        printf("with isotope envelope fitting, which is not fully supported yet. Please consider\n");
        printf("using a fixed number of isotope envelope peaks instead.\n");
    }
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
        else if (lk_Char.isLetter() || lk_Char == '*' || lk_Char == '^')
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


tk_IsotopeEnvelope k_QuantifierBase::heavyEnvelopeForPeptide(QString as_Peptide, tk_StringIntHash ak_StarAminoAcids)
{
    tk_IsotopeEnvelope lk_Result;

    foreach (QString ls_AminoAcid, mk_StarAminoAcids)
        if (!ak_StarAminoAcids.contains(ls_AminoAcid))
            ak_StarAminoAcids[ls_AminoAcid] = 0;
        
    tk_StringIntHash lk_AvailableStarAminoAcids;
    foreach (QString ls_AminoAcid, ak_StarAminoAcids.keys())
        lk_AvailableStarAminoAcids[ls_AminoAcid] = ak_StarAminoAcids[ls_AminoAcid];
    
    bool lb_First = true;
    
    for (int i = 0; i < as_Peptide.length(); ++i)
    {
        QString ls_AminoAcid = as_Peptide.mid(i, 1);
        tk_IsotopeEnvelope lk_UseEnvelope = mk_LightIsotopeEnvelopeForAminoAcid[ls_AminoAcid];
        if (mk_HeavyIsotopeEnvelopeForAminoAcid.contains(ls_AminoAcid))
        {
            if (mk_StarAminoAcids.contains(ls_AminoAcid))
            {
                if (lk_AvailableStarAminoAcids[ls_AminoAcid] > 0)
                {
                    lk_UseEnvelope = mk_HeavyIsotopeEnvelopeForAminoAcid[ls_AminoAcid];
                    --lk_AvailableStarAminoAcids[ls_AminoAcid];
                }
            }
            else
                lk_UseEnvelope = mk_HeavyIsotopeEnvelopeForAminoAcid[ls_AminoAcid];
        }
        if (lb_First)
            lk_Result = lk_UseEnvelope;
        else
            lk_Result = mk_IsotopeEnvelope.add(lk_Result, lk_UseEnvelope);
        lb_First = false;
    }
    return lk_Result;
    // :TODO: fix this for star amino acids
/*    QHash<QString, int> lk_Composition = compositionForPeptide(as_Peptide);
    
    QHash<QString, QHash<QString, int> > lk_CompositionForAminoAcid;
    for (int i = 0; i < as_Peptide.length(); ++i)
    {
        QString ls_AminoAcid = as_Peptide.mid(i, 1);
        QString ls_Id = QString();
        //if (mk_HeavyIsotopeEnvelopeForAminoAcid.contains(ls_AminoAcid) || mk_StarAminoAcids.contains(ls_AminoAcid))
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

    tk_StringIntHash lk_UsedStarAminoAcids;

    foreach (QString ls_AminoAcid, mk_StarAminoAcids)
        if (!ak_StarAminoAcids.contains(ls_AminoAcid))
            ak_StarAminoAcids[ls_AminoAcid] = 0;
        
    foreach (QString ls_AminoAcid, ak_StarAminoAcids.keys())
        lk_UsedStarAminoAcids[ls_AminoAcid] = 0;
    
    foreach (QString ls_Id, lk_CompositionForAminoAcid.keys())
    {
        k_IsotopeEnvelope* lk_UseEnvelope_ = &mk_IsotopeEnvelope;
        if (!ls_Id.isEmpty())
        {
            if (mk_StarAminoAcids.contains(ls_Id))
            {
                if (lk_UsedStarAminoAcids[ls_Id] < ak_StarAminoAcids[ls_Id])
                {
                    lk_UseEnvelope_ = &mk_HeavyIsotopeEnvelopeForAminoAcid[ls_Id];
                    ++lk_UsedStarAminoAcids[ls_Id];
                }
            }
            else
                lk_UseEnvelope_ = &mk_HeavyIsotopeEnvelopeForAminoAcid[ls_Id];
        }
        if (lb_First)
            lr_HeavyMass = lk_UseEnvelope_->isotopeEnvelopeForComposition(lk_CompositionForAminoAcid[ls_Id]);
        else
            lr_HeavyMass = mk_IsotopeEnvelope.add(lr_HeavyMass, lk_UseEnvelope_->isotopeEnvelopeForComposition(lk_CompositionForAminoAcid[ls_Id]));
        lb_First = false;
    }
    
    return lr_HeavyMass;*/
}


double k_QuantifierBase::lightMassForPeptide(QString as_Peptide)
{
    return mk_IsotopeEnvelope.massForComposition(compositionForPeptide(as_Peptide));
}


QString k_QuantifierBase::heavyEnvelopeTitle(tk_StringIntHash ak_StarAminoAcids)
{
    QString ls_Result;
    
    // lk_AminoAcidOrder contains all labeled amino acids, unconditional first, followed by star amino acids
    QStringList lk_AminoAcidOrder;
    foreach (QString ls_AminoAcid, mk_AminoAcidForDescription.values())
    {
        if ((!ak_StarAminoAcids.contains(ls_AminoAcid)) && (!lk_AminoAcidOrder.contains(ls_AminoAcid)))
            lk_AminoAcidOrder << ls_AminoAcid;
    }
    foreach (QString ls_AminoAcid, mk_AminoAcidForDescription.values())
    {
        if ((ak_StarAminoAcids.contains(ls_AminoAcid)) && (!lk_AminoAcidOrder.contains(ls_AminoAcid)))
            lk_AminoAcidOrder << ls_AminoAcid;
    }
    
    foreach (QString ls_Description, mk_AminoAcidForDescription.uniqueKeys())
    {
        if (!ls_Result.isEmpty())
            ls_Result += ", ";
        QString ls_AminoAcidScope;
        foreach (QString ls_AminoAcid, lk_AminoAcidOrder)
        {
            if (mk_AminoAcidForDescription.values(ls_Description).contains(ls_AminoAcid))
            {
                if (ak_StarAminoAcids.contains(ls_AminoAcid))
                {
                    if (ak_StarAminoAcids[ls_AminoAcid] == 0)
                        ls_AminoAcid = "";
                    else
                        ls_AminoAcid += QString("%1").arg(ak_StarAminoAcids[ls_AminoAcid]);
                }
                ls_AminoAcidScope += ls_AminoAcid;
            }
        }
        if (ls_AminoAcidScope.length() != 20)
            ls_Result += ls_AminoAcidScope + " ";
        ls_Result += ls_Description;
    }
    
    return ls_Result;
}


double k_QuantifierBase::nominalMassShiftForPeptide(QString as_Peptide, tk_StringIntHash ak_StarAminoAcids)
{
    QHash<QString, int> lk_Composition = compositionForPeptide(as_Peptide);
    
    double ld_MassShift = 0.0;
    
    tk_StringIntHash lk_UsedStarAminoAcids;
    
    foreach (QString ls_AminoAcid, mk_StarAminoAcids)
        if (!ak_StarAminoAcids.contains(ls_AminoAcid))
            ak_StarAminoAcids[ls_AminoAcid] = 0;
        
    foreach (QString ls_AminoAcid, ak_StarAminoAcids.keys())
        lk_UsedStarAminoAcids[ls_AminoAcid] = 0;
    
    QHash<QString, QHash<QString, int> > lk_CompositionForAminoAcid;
    for (int i = 0; i < as_Peptide.length(); ++i)
    {
        QString ls_AminoAcid = as_Peptide.mid(i, 1);
        double ld_Delta = 0.0;
        if (!mk_StarAminoAcids.contains(ls_AminoAcid))
            ld_Delta = mk_NominalMassShiftForAminoAcid[ls_AminoAcid];
        else
        {
            if (lk_UsedStarAminoAcids[ls_AminoAcid] < ak_StarAminoAcids[ls_AminoAcid])
            {
                ld_Delta = mk_NominalMassShiftForAminoAcid[ls_AminoAcid];
                ++lk_UsedStarAminoAcids[ls_AminoAcid];
            }
        }
        ld_MassShift += ld_Delta;
    }
    
    return ld_MassShift;
}
