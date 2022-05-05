#!/usr/bin/env python


import os
import time
import datetime
from urllib.error import HTTPError

import click
import google.auth
import googleapiclient.discovery

# without this it won't work! Have your own credentials.json file
os.environ['GOOGLE_APPLICATION_CREDENTIALS'] = '/home/opich/keys.json'


# [START create_instance]
def create_instance(compute, project, zone, name, bucket=None, preemptible=False,
                    tumoral=None, normal=None, cores=None):

    # FIXME
    source_disk_image = f"projects/{project}/global/images/image-4"

    # Configure the machine
    # FIXME
    machine_type = f"projects/{project}/zones/{zone}/machineTypes/custom-32-61440"

    # startup script
    startup_script = open(
        os.path.join(
            os.path.dirname(__file__), 'startup-script.sh'), 'r').read()

    config = {
        'name': name,
        'machineType': machine_type,

        # Specify the boot disk and the image to use as a source.
        'disks': [
            {
                'boot': True,
                'autoDelete': True,
                'initializeParams': {
                    'sourceImage': source_disk_image,
                    "diskSizeGb": "300"
                }
            }
        ],
        'enable-guest-attributes': True,
        # Specify a network interface with NAT to access the public
        # internet.
        'networkInterfaces': [{
            'network': 'global/networks/default',
            'accessConfigs': [
                {'type': 'ONE_TO_ONE_NAT', 'name': 'External NAT'}
            ]
        }],

        # Allow the instance to access cloud storage and logging.
        'serviceAccounts': [{
            'email': 'WRITEHEREYOUREMAILACCOUNT',
            'scopes': [
                "https://www.googleapis.com/auth/cloud-platform"
            ]
        }],

        # Metadata is readable from the instance and allows you to
        # pass configuration from deployment scripts to instances.
        'metadata': {
            'items': [{
                    # Startup script is automatically executed by the
                    # instance upon startup.
                    'key': 'startup-script',
                    'value': startup_script
                }, {
                    'key': 'tumoral',
                    'value': tumoral
                }, {
                    'key': 'normal',
                    'value': normal
                }, {
                    'key': 'bucket',
                    'value': bucket
                },
                {
                    'key': 'name',
                    'value': name
                },
                {
                    'key': 'zone',
                    'value': zone
                },
                {
                    'key': 'project',
                    'value': project
                },
                {
                    'key': 'cpus',
                    'value': cores
                },
                {
                    'key': 'enable-guest-attributes',
                    'value': True
                }]
        }
    }

    if preemptible:
        config['scheduling'] = [{'preemptible': True}]

    return compute.instances().insert(
        project=project,
        zone=zone,
        body=config).execute()


def delete_instance(compute, project, zone, name):
    return compute.instances().delete(
        project=project,
        zone=zone,
        instance=name).execute()


def wait_for_operation(compute, project, zone, operation):
    print('Waiting for operation to finish...')
    while True:
        result = compute.zoneOperations().get(
            project=project,
            zone=zone,
            operation=operation).execute()

        if result['status'] == 'DONE':
            print("done.")
            if 'error' in result:
                raise Exception(result['error'])
            return result

        time.sleep(1)


def delete_instance_full(compute, project, zone, instance_name):
    operation = delete_instance(compute, project, zone, instance_name)
    wait_for_operation(compute, project, zone, operation['name'])


def main(project=None, zone=None, instance_name=None, **kwargs):
    start_time = datetime.datetime.now()
    compute = googleapiclient.discovery.build('compute', 'v1')

    # credentials, project_ = google.auth.default()

    print('Creating instance.')

    operation = create_instance(compute, project, zone, instance_name, **kwargs)
    wait_for_operation(compute, project, zone, operation['name'])

    print('Instance created!!')
    # if it is preempted it will
    running = True
    error = True
    old_status_pipeline = 'TEST'
    while running is True:
        # if the instance is dead
        try:
            metadata = compute.instances().get(
                instance=instance_name,
                zone=zone, project=project
            ).execute(num_retries=10)

        except HTTPError:
            print('Deleting instance because VM is preempted')
            running = False
        else:

            time.sleep(10)
            for elements in metadata['metadata']['items']:
                if elements['key'] == 'STATUSPIPELINE':
                    current_status_pipeline = elements['value']
                    if current_status_pipeline != old_status_pipeline:
                        new_time = datetime.datetime.now() - start_time
                        print(new_time, current_status_pipeline)
                        old_status_pipeline = current_status_pipeline

            if metadata['status'] == 'RUNNING':
                for elements in metadata['metadata']['items']:
                    if elements['key'] == 'finished_task':
                        if elements['value'] == 'succeed':
                            running = False
                            error = False
                        else:
                            print("Something has failed...")
                            running = False
            else:
                print(" the instance is not running")
                running = False

    print('Deleting instance because VM job has finished')
    delete_instance_full(compute, project, zone, instance_name)

    if error:
        quit(-1)


@click.command()
@click.option('-n', '--normal', required=True, help='Path to normal bam')
@click.option('-t', '--tumoral', required=True, help='Path to tumoral bam')
@click.option('-b', '--bucket', required=True, help='Bucket for the output')
@click.option('-c', '--cores', default=1, help='Cores to use')
@click.option('-p', '--project', help='Project name')
@click.option('-z', '--zone', default='europe-west4-b', help='Compute Engine zone to deploy to.')
@click.option('-i', '--name', 'instance_name', default='demo-instance', help='Instance name')
@click.option('--preemptible', is_flag=True)
def cli(**kwargs):
    main(**kwargs)


if __name__ == '__main__':
    cli()
