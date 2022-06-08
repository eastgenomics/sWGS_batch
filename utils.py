import datetime

import dxpy


def get_run_datetime():
    """ Return datetime as string in the following format: YYMMDD-HHmm

    Returns:
        str: Datetime in string format
    """

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


def create_dnanexus_links(list_dnanexus_ids, app_handler, input_name):
    """ Return dnanexus links when passed a iterable of dnanexus ids

    Args:
        list_dnanexus_ids (list): List of dnanexus ids

    Returns:
        list: List of dnanexus links
    """

    list_dnanexus_links = []

    for dnanexus_id in list_dnanexus_ids:
        dnanexus_link = {"$dnanexus_link": dnanexus_id}
        list_dnanexus_links.append(dnanexus_link)

    input_type = get_input_type(app_handler, input_name)

    if input_type == "array:file":
        return list_dnanexus_links
    elif input_type == "file":
        if len(list_dnanexus_links) == 1:
            return list_dnanexus_links[0]
        else:
            raise Exception((
                f"Input type for {input_name} is 'file' but multiple inputs "
                f"have been gathered: {list_dnanexus_ids}"
            ))


def get_input_type(app_handler, input_name):
    """ Get the input type given the input name

    Args:
        app_handler (DXApp): DXApp for WisecondorX
        input_name (str): Input name

    Returns:
        str: Type of the input
    """

    stage_input = None

    for stage_input in app_handler.describe()["inputSpec"]:
        if stage_input["name"] == input_name:
            return stage_input["class"]

    return stage_input
