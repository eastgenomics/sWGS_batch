import datetime

import dxpy


def get_run_datetime():
    now = datetime.datetime.now()
    run_datetime = now.strftime("%y%m%d-%H%M")
    return run_datetime


def get_workflow_stage_info(workflow_id):
    """ Get the workflow stage info i.e. stage id, app id and app name
    Args:
        workflow_id (str): Workflow id
    Returns:
        dict: Dict of stage id to app info
    """

    workflow_description_json = dxpy.describe(workflow_id)

    stages = {}

    # go through the workflow json description and select stage,
    # app_id and app_name
    for stage in workflow_description_json['stages']:
        # gather app id and app name of the stage
        app_id = stage['executable']
        app_name = dxpy.describe(app_id)['name']
        stages[stage['id']] = {
            "app_id": app_id, "app_name": app_name
        }

    return stages


def get_stage_output_folders(workflow_stage_info):
    """ Add specific app output directory to final command line
    Args:
        workflow_stage_info (dict): Dict containing the stage to app info
    Returns:
        str: String containing the DNAnexus option to specify apps output directory
    """

    stage_folders = {}

    for stage_id, stage_info in sorted(workflow_stage_info.items()):
        stage_folders[stage_id] = stage_info["app_name"]
    return stage_folders


def create_dnanexus_links(list_dnanexus_ids):
    list_dnanexus_links = []

    for dnanexus_id in list_dnanexus_ids:
        dnanexus_link = {"$dnanexus_link": dnanexus_id}
        list_dnanexus_links.append(dnanexus_link)

    if len(list_dnanexus_links) == 1:
        return list_dnanexus_links[0]
    else:
        return list_dnanexus_links
