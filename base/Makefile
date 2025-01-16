# multipass convenience targets (host side)
#  -- multipass is a nice and easy ubuntu vm provider for macs
#  -- mutlipass enables singularity development and testing on macs,
#  -- that is otherwise impossible directly on the ARM architecture

.PHONY: mp_launch mp_start mp_stop mp_rm mp_shell mp_transfer

MP_INSTANCE_NAME := singer

DOCKER_IMAGE_HOST := /Users/pasha/var/docker_images/base.tar
LICENSE_HOST := /Users/pasha/licenses/chemaxon/license.cxl

DOCKER_IMAGE_MP = /home/ubuntu/base.docker.img.tar
LICENSE_MP = /home/ubuntu/license.cxl


mp_launch: vm_spec.yaml
	multipass launch \
		--cloud-init $< \
		-d 8G -m 4G -c 2 \
		-n $(MP_INSTANCE_NAME) \
		noble

mp_transfer:
	multipass transfer -p $(DOCKER_IMAGE_HOST) $(MP_INSTANCE_NAME):$(DOCKER_IMAGE_MP)
	multipass transfer -p $(LICENSE_HOST) $(MP_INSTANCE_NAME):$(LICENSE_MP)
	multipass transfer -p Makefile $(MP_INSTANCE_NAME):/home/ubuntu/Makefile

mp_shell: mp_transfer
	multipass shell $(MP_INSTANCE_NAME)

mp_start:
	multipass start $(MP_INSTANCE_NAME)

mp_stop:
	multipass stop $(MP_INSTANCE_NAME)

mp_rm: mp_stop
	multipass delete -p $(MP_INSTANCE_NAME)

# docker convenience targets (host side)

.PHONY: docker_build docker_push docker_run


DOCKER_IMAGE_NAME := synbiorox/base
DOCKER_IMAGE_VERSION := 0.1.0
DOCKER_IMAGE_TAG = $(DOCKER_IMAGE_NAME):$(DOCKER_IMAGE_VERSION)

docker_build:
	docker build . -t $(DOCKER_IMAGE_TAG)


docker_push:
	docker push $(DOCKER_IMAGE_TAG)


$(HOST_DOCKER_IMAGES)/%.tar:
	docker save $(DOCKER_IMAGE_TAG) -o $@

docker_save: $(HOST_DOCKER_IMAGES)/base.tar

docker_run:
	docker run \
		--rm \
		--mount type=bind,source=$(HOST_LICENSES),target=/opt/chemaxon/jchemsuite/licenses \
		-it $(DOCKER_IMAGE_TAG) \
		bash

# singularity convenience (the main thing we are trying to develop and test)
#   -- meant to be used from inside the multipass vm or from quest

.PHONY: singularity_build singularity_test_cxcalc

IMAGE_SINGULARITY = /home/ubuntu/base.sif
LICENSE_SINGULARITY = /home/ubuntu/license.cxl

%.sif:
	singularity build $@ docker-archive://$(DOCKER_IMAGE_MP)

singularity_build: $(IMAGE_SINGULARITY)

singularity_test_cxcalc: $(IMAGE_SINGULARITY)
	singularity run \
		--env CHEMAXON_LICENSE_URL=$(LICENSE_SINGULARITY) \
		--bind $(LICENSE_MP):$(LICENSE_SINGULARITY):ro \
		--bind /home/ubuntu/audit:/opt/chemaxon/jchemsuite/audit:rw \
		$< \
		cxcalc pKa "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
