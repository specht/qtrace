#include "Tango.h"


QString tangoMix(const QString& a, const QString& b, unsigned char auc_Factor)
{
    bool lb_Ok;
    int ar = a.mid(1, 2).toInt(&lb_Ok, 16);
    int ag = a.mid(3, 2).toInt(&lb_Ok, 16);
    int ab = a.mid(5, 2).toInt(&lb_Ok, 16);
    int br = b.mid(1, 2).toInt(&lb_Ok, 16);
    int bg = b.mid(3, 2).toInt(&lb_Ok, 16);
    int bb = b.mid(5, 2).toInt(&lb_Ok, 16);
    int mr = (ar * auc_Factor + br * (255 - auc_Factor)) / 255;
    int mg = (ag * auc_Factor + bg * (255 - auc_Factor)) / 255;
    int mb = (ab * auc_Factor + bb * (255 - auc_Factor)) / 255;
    return QString("#%1%2%3")
        .arg(mr, 2, 16, QChar('0'))
        .arg(mg, 2, 16, QChar('0'))
        .arg(mb, 2, 16, QChar('0'));
}
