# multipass convenience targets (host side)
#  -- multipass is a nice and easy ubuntu vm provider for macs
#  -- mutlipass enables singularity development and testing on macs,
#  -- that is otherwise impossible directly on the ARM architecture

.PHONY: mp_launch mp_start mp_stop mp_rm mp_shell mp_transfer

MP_INSTANCE_NAME := singer

DOCKER_IMAGE_HOST := $(HOME)/var/docker_images/base.tar
LICENSE_HOST := $(HOME)/licenses/chemaxon/license.cxl

DOCKER_IMAGE_MP = $(HOME)/base.docker.img.tar
LICENSE_MP = $(HOME)/license.cxl


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

platform_amd64 = --platform=linux/amd64

docker_build.amd64:
	docker build . -t $(DOCKER_IMAGE_TAG) $(platform_amd64)


docker_push.amd64:
	docker push $(DOCKER_IMAGE_TAG) $(platform_amd64)


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

.PHONY: singularity_build singularity_pull singularity_test_cxcalc

IMAGE_SINGULARITY = $(HOME)/bottle_0.1.0.sif

license_singularity = /home/bottler/license.cxl
audit_dir_host = $(HOME)/var/chemaxon/audit

%.sif:
	singularity build $@ docker-archive://$(DOCKER_IMAGE_MP)

singularity_build: $(IMAGE_SINGULARITY)

# -- quest main

singularity_pull:
	singularity pull docker://$(DOCKER_IMAGE_TAG)

$(audit_dir_host):
	mkdir -p $@

singularity_test_cxcalc: $(IMAGE_SINGULARITY) $(audit_dir_host)
	singularity run \
		--env CHEMAXON_LICENSE_URL=$(license_singularity) \
		--bind $(LICENSE_MP):$(license_singularity):ro \
		--bind $(audit_dir_host):/opt/chemaxon/jchemsuite/audit:rw \
		$< \
		cxcalc pKa "Cn1cnc2c1c(=O)n(C)c(=O)n2C"
