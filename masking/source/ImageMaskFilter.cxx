#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkConfigure.h"
#include "itkMaskImageFilter.h"
#include "itkImageRegionIterator.h"


typedef itk::Image<float, 3>  ImageType;
typedef itk::Image<float, 3>  OutImageType;
typedef itk::ImageFileReader< ImageType > ReaderType;
typedef itk::ImageFileWriter< OutImageType > WriterType;


int main(int argc, char * argv[])
{
	ReaderType::Pointer reader1 = ReaderType::New();//ChestFixed
	ReaderType::Pointer reader2 = ReaderType::New();//LungFixed
	ReaderType::Pointer reader3 = ReaderType::New();//ChestMoving
	ReaderType::Pointer reader4 = ReaderType::New();//LungMoving

	WriterType::Pointer writer1 = WriterType::New();//Chestfixed+LungFixed=MaskedFixed
	WriterType::Pointer writer2= WriterType::New(); //ChestMoving+LungMoving =MaskedMoving

	//Parameters
	reader1->SetFileName(argv[1]);
	reader2->SetFileName(argv[2]);
	reader3->SetFileName(argv[3]);
	reader4->SetFileName(argv[4]);
	writer1->SetFileName(argv[5]);
	writer2->SetFileName(argv[6]);
	
	//Pipeline
	try
	{
		reader1->Update();
		reader2->Update();
		reader3->Update();
		reader4->Update();
	}
	catch (itk::ExceptionObject &err)
	{
		std::cout << "Problems reading input image" << std::endl;
		std::cerr << "ExceptionObject caught!" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}


	typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter1 = MaskFilterType::New();
	maskFilter1->SetInput(reader1->GetOutput());
	maskFilter1->SetMaskImage(reader2->GetOutput());
	
	writer1->SetInput(maskFilter1->GetOutput());
	
	
	typedef itk::MaskImageFilter< ImageType, ImageType > MaskFilterType;
	MaskFilterType::Pointer maskFilter2 = MaskFilterType::New();
	maskFilter2->SetInput(reader3->GetOutput());
	maskFilter2->SetMaskImage(reader4->GetOutput());

	writer2->SetInput(maskFilter2->GetOutput());

	try
	{
		writer1->Update();
		writer2->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cout << "ExceptionObject caught !" << std::endl;
		std::cout << err << std::endl;
		return EXIT_FAILURE;
	}


	return EXIT_SUCCESS;


}