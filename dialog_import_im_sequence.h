#ifndef DIALOG_IMPORT_IM_SEQUENCE_H
#define DIALOG_IMPORT_IM_SEQUENCE_H
#include "ui_dialog_import_im_sequence.h"

#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"

struct ImportImgPara
{
    V3DLONG startImg;
    V3DLONG endImg;
    V3DLONG countImg;
    QString filt;
    V3DLONG inc;
    int packType;
};

class dialogImportImSequence : public QDialog, public Ui_import_im_seq
{
    Q_OBJECT

public:
    dialogImportImSequence()
    {
        setupUi(this);
        connect(OK_button, SIGNAL(clicked()), this, SLOT(accept()));
        connect(Cancel_button, SIGNAL(clicked()), this, SLOT(reject()));

        //initialization codes for read directory
        num_im_SpBox->setMaximum(10000000);
        start_id_SpBox->setMaximum(10000000);
        start_id_SpBox->setMinimum(1);
        increment_SpBox->setMaximum(10000000);
        increment_SpBox->setMinimum(1);
        end_id_SpBox->setMaximum(10000000);
        end_id_SpBox->setMinimum(1);

        start_id_SpBox->setValue(1);
        increment_SpBox->setValue(1);
        end_id_SpBox->setValue(1);
        LE_name_key_word->setText("");

        //comboPack->clear();
        //comboPack->insertItem(0, "Pack images in 'Z' dimension");
        //comboPack->insertItem(1, "Pack images in 'Channel' dimension");
        //comboPack->setCurrentIndex(0);
    }

    void fetchData(ImportImgPara *p)
    {
        if (!p) return;
        p->countImg = num_im_SpBox->value();
        p->startImg = start_id_SpBox->value();
        p->endImg = end_id_SpBox->value();
        p->inc = increment_SpBox->value();
        p->filt = LE_name_key_word->text();
        //this equals to setting p->packType is set to one, if using combo list.
        p->packType = TIME_PACK_C; //for 3d time series data, we all concatenate them in Channel dim
    }

private:
    QStringList mystringlist;

};

#endif // DIALOG_IMPORT_IM_SEQUENCE_H
