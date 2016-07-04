#include "itkImageRegistrationMethod.h"
#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkTimeProbesCollectorBase.h"
#include "itkMemoryProbesCollectorBase.h"
#include "itkBSplineTransform.h"
#include "itkRegularStepGradientDescentOptimizer.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSquaredDifferenceImageFilter.h"
#include "itkTransformFileReader.h"

int main(int argc, char *argv[])
{
	if (argc < 4

	{
		std::cerr << "Missing Parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " fixedImageFile  movingImageFile outputImagefile  ";
		std::cerr << " [differenceOutputfile] [differenceBeforeRegistration] ";
		std::cerr << " [deformationField] ";
		std::cerr << " [useExplicitPDFderivatives ] [useCachingBSplineWeights ] ";
		std::cerr << " [filenameForFinalTransformParameters] ";
		std::cerr << " [maximumStepLength] [maximumNumberOfIterations]";
		std::cerr << std::endl;
		return EXIT_FAILURE;
	}

	const    unsigned int    ImageDimension = 3;
	typedef  signed short    PixelType;

	typedef itk::Image< PixelType, ImageDimension >  FixedImageType;
	typedef itk::Image< PixelType, ImageDimension >  MovingImageType;

	const unsigned int SpaceDimension = ImageDimension;
	const unsigned int SplineOrder = 3;
	typedef double CoordinateRepType;

	typedef itk::BSplineTransform<
		CoordinateRepType,
		SpaceDimension,
		SplineOrder >     TransformType;

	typedef itk::RegularStepGradientDescentOptimizer       OptimizerType;

	typedef itk::MattesMutualInformationImageToImageMetric<
		FixedImageType,
		MovingImageType >    MetricType;

	typedef itk::LinearInterpolateImageFunction<
		MovingImageType,
		double          >    InterpolatorType;

	typedef itk::ImageRegistrationMethod<
		FixedImageType,
		MovingImageType >    RegistrationType;

	MetricType::Pointer         metric = MetricType::New();
	OptimizerType::Pointer      optimizer = OptimizerType::New();
	InterpolatorType::Pointer   interpolator = InterpolatorType::New();
	RegistrationType::Pointer   registration = RegistrationType::New();

	registration->SetMetric(metric);
	registration->SetOptimizer(optimizer);
	registration->SetInterpolator(interpolator);

	TransformType::Pointer  transform = TransformType::New();
	registration->SetTransform(transform);

	typedef itk::ImageFileReader< FixedImageType  > FixedImageReaderType;
	typedef itk::ImageFileReader< MovingImageType > MovingImageReaderType;

	FixedImageReaderType::Pointer  fixedImageReader = FixedImageReaderType::New();
	MovingImageReaderType::Pointer movingImageReader = MovingImageReaderType::New();

	fixedImageReader->SetFileName(argv[1]);
	movingImageReader->SetFileName(argv[2]);

	FixedImageType::ConstPointer fixedImage = fixedImageReader->GetOutput();

	registration->SetFixedImage(fixedImage);
	registration->SetMovingImage(movingImageReader->GetOutput());

	fixedImageReader->Update();

	FixedImageType::RegionType fixedRegion = fixedImage->GetBufferedRegion();

	registration->SetFixedImageRegion(fixedRegion);

	unsigned int numberOfGridNodesInOneDimension = 5;



	TransformType::PhysicalDimensionsType   fixedPhysicalDimensions;
	TransformType::MeshSizeType             meshSize;
	TransformType::OriginType               fixedOrigin;

	for (unsigned int i = 0; i< SpaceDimension; i++)
	{
		fixedOrigin[i] = fixedImage->GetOrigin()[i];
		fixedPhysicalDimensions[i] = fixedImage->GetSpacing()[i] *
			static_cast<double>(
			fixedImage->GetLargestPossibleRegion().GetSize()[i] - 1);
	}
	meshSize.Fill(numberOfGridNodesInOneDimension - SplineOrder);

	transform->SetTransformDomainOrigin(fixedOrigin);
	transform->SetTransformDomainPhysicalDimensions(
		fixedPhysicalDimensions);
	transform->SetTransformDomainMeshSize(meshSize);
	transform->SetTransformDomainDirection(fixedImage->GetDirection());

	typedef TransformType::ParametersType     ParametersType;

	const unsigned int numberOfParameters =
		transform->GetNumberOfParameters();

	ParametersType parameters(numberOfParameters);

	parameters.Fill(0.0);

	transform->SetParameters(parameters);

	registration->SetInitialTransformParameters(transform->GetParameters());
	
	parameters.Fill(0.0);

	transform->SetParameters(parameters);

	registration->SetInitialTransformParameters(transform->GetParameters());




	//  Next we set the parameters of the RegularStepGradientDescentOptimizer object.
	
	optimizer->SetMaximumStepLength(15.0); // Default =10 , the maximum value is 20, with increase in the maximum step length the accuracy and computation time is found to decrease.
	optimizer->SetMinimumStepLength(0.009);//Default value = 0.01, Computation time increases beyond 0.009 value. Better convergence is observed at this point.

	optimizer->SetRelaxationFactor(0.65);// default= 0.7, since the given image is bit less noisy. It was decreased for 0.65. Decreasing further would decrease the step size. 
	optimizer->SetNumberOfIterations(85);// Default value is 50, runs out of memory for > 90


	// Optionally, get the step length from the command line arguments
	if (argc > 12)
	{
		optimizer->SetMaximumStepLength(atof(argv[12]));
	}
	// Optionally, get the number of iterations from the command line arguments
	if (argc > 13)
	{
		optimizer->SetNumberOfIterations(atoi(argv[13]));
	}

	// Create the Command observer and register it with the optimizer.
	CommandIterationUpdate::Pointer observer = CommandIterationUpdate::New();
	optimizer->AddObserver(itk::IterationEvent(), observer);

	metric->SetNumberOfHistogramBins(50);

	const unsigned int numberOfSamples =
		static_cast<unsigned int>(fixedRegion.GetNumberOfPixels() * 20.0 / 100.0);

	metric->SetNumberOfSpatialSamples(numberOfSamples);
	metric->ReinitializeSeed(76926294);

	if (argc > 7)
	{
		// Define whether to calculate the metric derivative by explicitly
		// computing the derivatives of the joint PDF with respect to the Transform
		// parameters, or doing it by progressively accumulating contributions from
		// each bin in the joint PDF.
		metric->SetUseExplicitPDFDerivatives(atoi(argv[7]));
	}

	if (argc > 8)
	{
		// Define whether to cache the BSpline weights and indexes corresponding to
		// each one of the samples used to compute the metric. Enabling caching will
		// make the algorithm run faster but it will have a cost on the amount of memory
		// that needs to be allocated. This option is only relevant when using the
		// BSplineTransform.
		metric->SetUseCachingOfBSplineWeights(atoi(argv[8]));
	}

	// Add time and memory probes
	itk::TimeProbesCollectorBase chronometer;
	itk::MemoryProbesCollectorBase memorymeter;

	std::cout << std::endl << "Starting Registration" << std::endl;

	try
	{
		memorymeter.Start("Registration");
		chronometer.Start("Registration");

		registration->Update();

		chronometer.Stop("Registration");
		memorymeter.Stop("Registration");

		std::cout << "Optimizer stop condition = "
			<< registration->GetOptimizer()->GetStopConditionDescription()
			<< std::endl;
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	OptimizerType::ParametersType finalParameters =
		registration->GetLastTransformParameters();
	transform->SetParameters(finalParameters);

	typedef itk::ResampleImageFilter<
		MovingImageType,
		FixedImageType >    ResampleFilterType;

	ResampleFilterType::Pointer resample = ResampleFilterType::New();

	resample->SetTransform(transform);
	resample->SetInput(movingImageReader->GetOutput());

	resample->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
	resample->SetOutputOrigin(fixedImage->GetOrigin());
	resample->SetOutputSpacing(fixedImage->GetSpacing());
	resample->SetOutputDirection(fixedImage->GetDirection());

	// This value is set to zero in order to make easier to perform
	// regression testing in this example. However, for didactic
	// exercise it will be better to set it to a medium gray value
	// such as 100 or 128.
	resample->SetDefaultPixelValue(0);

	typedef  signed short  OutputPixelType;

	typedef itk::Image< OutputPixelType, ImageDimension > OutputImageType;

	typedef itk::CastImageFilter<
		FixedImageType,
		OutputImageType > CastFilterType;

	typedef itk::ImageFileWriter< OutputImageType >  WriterType;

	WriterType::Pointer      writer = WriterType::New();
	CastFilterType::Pointer  caster = CastFilterType::New();

	writer->SetFileName(argv[3]);

	caster->SetInput(resample->GetOutput());
	writer->SetInput(caster->GetOutput());

	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & err)
	{
		std::cerr << "ExceptionObject caught !" << std::endl;
		std::cerr << err << std::endl;
		return EXIT_FAILURE;
	}

	typedef itk::SquaredDifferenceImageFilter<
		FixedImageType,
		FixedImageType,
		OutputImageType > DifferenceFilterType;

	DifferenceFilterType::Pointer difference = DifferenceFilterType::New();

	WriterType::Pointer writer2 = WriterType::New();
	writer2->SetInput(difference->GetOutput());

	// Compute the difference image between the
	// fixed and resampled moving image.
	if (argc > 4)
	{
		difference->SetInput1(fixedImageReader->GetOutput());
		difference->SetInput2(resample->GetOutput());
		writer2->SetFileName(argv[4]);
		try
		{
			writer2->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}

	// Compute the difference image between the
	// fixed and moving image before registration.
	if (argc > 5)
	{
		writer2->SetFileName(argv[5]);
		difference->SetInput1(fixedImageReader->GetOutput());
		difference->SetInput2(movingImageReader->GetOutput());
		try
		{
			writer2->Update();
		}
		catch (itk::ExceptionObject & err)
		{
			std::cerr << "ExceptionObject caught !" << std::endl;
			std::cerr << err << std::endl;
			return EXIT_FAILURE;
		}
	}

	// Generate the explicit deformation field resulting from
	// the registration.
	if (argc > 6)
	{

		typedef itk::Vector< float, ImageDimension >      VectorType;
		typedef itk::Image< VectorType, ImageDimension >  DisplacementFieldType;

		DisplacementFieldType::Pointer field = DisplacementFieldType::New();
		field->SetRegions(fixedRegion);
		field->SetOrigin(fixedImage->GetOrigin());
		field->SetSpacing(fixedImage->GetSpacing());
		field->SetDirection(fixedImage->GetDirection());
		field->Allocate();

		typedef itk::ImageRegionIterator< DisplacementFieldType > FieldIterator;
		FieldIterator fi(field, fixedRegion);

		fi.GoToBegin();

		TransformType::InputPointType  fixedPoint;
		TransformType::OutputPointType movingPoint;
		DisplacementFieldType::IndexType index;

		VectorType displacement;

		while (!fi.IsAtEnd())
		{
			index = fi.GetIndex();
			field->TransformIndexToPhysicalPoint(index, fixedPoint);
			movingPoint = transform->TransformPoint(fixedPoint);
			displacement = movingPoint - fixedPoint;
			fi.Set(displacement);
			++fi;
		}

		typedef itk::ImageFileWriter< DisplacementFieldType >  FieldWriterType;
		FieldWriterType::Pointer fieldWriter = FieldWriterType::New();

		fieldWriter->SetInput(field);

		fieldWriter->SetFileName(argv[6]);
		try
		{
			fieldWriter->Update();
		}
		catch (itk::ExceptionObject & excp)
		{
			std::cerr << "Exception thrown " << std::endl;
			std::cerr << excp << std::endl;
			return EXIT_FAILURE;
		}
	}

	// save the transform parameters in a file
	if (argc > 9)
	{
		std::ofstream parametersFile;
		parametersFile.open(argv[9]);
		parametersFile << finalParameters << std::endl;
		parametersFile.close();
	}

	return EXIT_SUCCESS;
}
