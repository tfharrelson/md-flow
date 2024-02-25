from setuptools import setup

package_data = {
    'md_flow': [
        'md_inputs/*.mdp'
    ]
}

setup(package_data=package_data)
