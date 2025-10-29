MAIN_DIR=$1;
GCP_BUCKET_ADDRESS=$2;

#---------------------------------------------------------------------------
# MOVE TO BUCKET 
#---------------------------------------------------------------------------

echo "INFO: Moving to GCP bucket at $GCP_BUCKET_ADDRESS ..."

# move directory into bucket
echo "Copying ${MAIN_DIR} into bucket ${GCP_BUCKET_ADDRESS}"
gsutil -m cp -r ${MAIN_DIR} ${GCP_BUCKET_ADDRESS}

STATUS=$?

if [ $STATUS -eq 0 ]; then
    echo "INFO: successfully copied resources to GCP bucket ${GCP_BUCKET_ADDRESS}, removing local copy ..."
    rm -r ${MAIN_DIR}
    echo "INFO: Done initializing resources, will be found at ${GCP_BUCKET_ADDRESS}"
else
    echo "ERROR: Failed to copy to bucket $GCP_BUCKET_ADDRESS"
    exit $STATUS
fi

