PLATFORM=${1:-"Unknown"};
BUCKET_NAME=${2:-"Unknown"};
NUM_THREADS=${3:-"10"};
RESOURCES_DIR=${4:-"Resources"};
STAR_INDEX_DIR=${5:-"STAR2.7.10b"};
SALMON_INDEX_DIR=${6:-"Salmon1.10.2"};
MAIN_DIR=${7:-"main"};


#---------------------------------------------------------------------------
# Assertions
#---------------------------------------------------------------------------
if [ $PLATFORM != 'GCP' ] && [ $PLATFORM != 'AWS' ]; then
    >&2 echo "Unknown platform $PLATFORM, choose GCP or AWS"
    exit 1
fi

if [ $BUCKET_NAME == 'Unknown' ]; then
    >&2 echo "Bucket name unset, please provide bucket name as second argument (not including s3:// or gs://)"
    exit 1
fi


#---------------------------------------------------------------------------
# Constants
#---------------------------------------------------------------------------
# get script dir
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"



#---------------------------------------------------------------------------
# DOWNLOAD RESOURCES
#---------------------------------------------------------------------------

# download resources
echo "INFO: Downloading resources ..."
sh ${SCRIPT_DIR}/resources/download_resources.sh $MAIN_DIR/$RESOURCES_DIR
echo "INFO: Downloading resources done"

#---------------------------------------------------------------------------
# ADD SITE FILES
#---------------------------------------------------------------------------
echo "INFO: Fixing known sites file, see doi.org/10.1038/s41467-022-28841-4 for sites information ..."
# fix Gabay et al sites list + copy Pinto et al sites list
bash ${SCRIPT_DIR}/resources/sites/copy_and_index_sites.bash $MAIN_DIR/$RESOURCES_DIR
echo "INFO: Fixing known sites file done"
echo "INFO: Copying INI file ..."
# copy INI file
sh ${SCRIPT_DIR}/resources/sites/copy_cmp_ini_files.sh $MAIN_DIR/$RESOURCES_DIR
echo "INFO: Copying INI file done"

#---------------------------------------------------------------------------
# ADD REGION FILES
#---------------------------------------------------------------------------
# copy and index inverted repeats in 3'UTR files
echo "INFO: Copying and indexing CEI regions files, see TODO DOI for more information ..."
bash ${SCRIPT_DIR}/resources/regions/copy_and_index_regions.bash $MAIN_DIR/$RESOURCES_DIR
echo "INFO: Copying and indexing CEI regions files done"

#---------------------------------------------------------------------------
# CHECK FOR MISSING FILES
#---------------------------------------------------------------------------
echo "INFO: Checking for missing files ..."
bash ${SCRIPT_DIR}/helpers/check_for_missing_files.sh

#---------------------------------------------------------------------------
# CREATE INDICES
#---------------------------------------------------------------------------
# index genomes for STAR
echo "INFO: Generating STAR indices ..."
sh ${SCRIPT_DIR}/resources/STAR-genome-generate.sh $MAIN_DIR/$RESOURCES_DIR $MAIN_DIR/$STAR_INDEX_DIR $NUM_THREADS
echo "INFO: Generating STAR indices done"
# index transcriptome for Salmon
echo "INFO: Generating Salmon index ..."
sh ${SCRIPT_DIR}/resources/Salmon-genome-generate.sh $MAIN_DIR/$RESOURCES_DIR $MAIN_DIR/$SALMON_INDEX_DIR $NUM_THREADS
echo "INFO: Generating Salmon index done"

#---------------------------------------------------------------------------
# MOVE TO BUCKET 
#---------------------------------------------------------------------------

echo "INFO: Moving to bucket for $PLATFORM ..."

if [ $PLATFORM = 'GCP' ]; then
    bash ${SCRIPT_DIR}/helpers/move_to_GCP_bucket.bash $MAIN_DIR "gs://${BUCKET_NAME}"
elif [ $PLATFORM = 'AWS' ]; then
    bash ${SCRIPT_DIR}/helpers/move_to_AWS_bucket.bash $MAIN_DIR "s3://${BUCKET_NAME}"
else
    >&2 echo "Unknown platform $PLATFORM"
    exit 1
fi

echo "INFO: Moving to bucket for $PLATFORM done"
echo "bye bye"


