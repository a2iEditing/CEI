MAIN_DIR=$1;
AWS_BUCKET_ADDRESS=$2;


#---------------------------------------------------------------------------
# MOVE TO BUCKET 
#---------------------------------------------------------------------------

echo "INFO: Moving to AWS bucket at $GCP_BUCKET_ADDRESS ..."

# move directory into bucket
echo "Copying ${MAIN_DIR} into bucket ${AWS_BUCKET_ADDRESS}"
# move directory into bucket
aws s3 mv ${MAIN_DIR} ${AWS_BUCKET_ADDRESS}/main --recursive

STATUS=$?

if [ $STATUS -eq 0 ]; then
    echo "INFO: successfully copied resources to AWS S3 bucket ${AWS_BUCKET_ADDRESS}, removing local copy ..."
    rm -r ${MAIN_DIR}
    echo "INFO: Done initializing resources, will be found at ${AWS_BUCKET_ADDRESS}"
else
    echo "ERROR: Failed to copy to bucket $AWS_BUCKET_ADDRESS"
    exit $STATUS
fi

