#ifndef DIALOG_IMPORT_IM_SEQUENCE_H
#define DIALOG_IMPORT_IM_SEQUENCE_H
//#include "ui_.h"

//#include "src_3rd/basic_c_fun/v3d_basicdatatype.h"

//struct ImportImgPara
//{
//    V3DLONG startImg;
//    V3DLONG endImg;
//    V3DLONG countImg;
//	QString filt;
//    V3DLONG inc;
//	int packType;
//};

//class import_images_tool_Dialog : public QDialog, public Ui_import_images_tool
//{
//    Q_OBJECT

//public:
//    import_images_tool_Dialog(QString folderName)
//	{
//		setupUi(this);

//		connect(okButton, SIGNAL(clicked()), this, SLOT(accept()));
//		connect(cancelButton, SIGNAL(clicked()), this, SLOT(reject()));

//		//initialization codes for read directory
//        numimgBox->setMaximum(10000000);
//        startimgBox->setMaximum(10000000);
//		startimgBox->setMinimum(1);
//        incBox->setMaximum(10000000);
//        incBox->setMinimum(1);
//        endimgBox->setMaximum(10000000);
//        endimgBox->setMinimum(1);

//		startimgBox->setValue(1);
//		incBox->setValue(1);
//        endimgBox->setValue(1);
//        filterEdit->setText("");

//        //comboPack->clear();
//        //comboPack->insertItem(0, "Pack images in 'Z' dimension");
//        //comboPack->insertItem(1, "Pack images in 'Channel' dimension");
//        //comboPack->setCurrentIndex(0);
//	}

//	void fetchData(ImportImgPara *p)
//	{
//		if (!p) return;
//		p->countImg = numimgBox->value();
//		p->startImg = startimgBox->value();
//        p->endImg = endimgBox->value();
//        p->inc = incBox->value();
//		p->filt = filterEdit->text();
//        //p->packType = comboPack->currentIndex();
//        p->packType = 1; //for 3d time series data, we all concatenate them in Channel dim
//	}

//private:
//	QStringList mystringlist;

//};

#endif // DIALOG_IMPORT_IM_SEQUENCE_H
