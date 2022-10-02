"""
Output autoupdate PR text
"""

import argparse
import yaml


parser = argparse.ArgumentParser()
parser.add_argument('--repo', required=True, help='Tool repo')
parser.add_argument('--old_version', required=True, help='Old version')
parser.add_argument('--new_version', required=True, help='New version')
parser.add_argument('--shed', required=True, help='Location of .shed.yml file input.')
parser.add_argument('--out', required=True, help='Output file.')
args = parser.parse_args()

update = f"from version {args.old_version} to {args.new_version}"

text = []
text.append(f"Hello! This is an automated update of the following tool: **{args.repo}**. I created this PR because I think the tool's main dependency is out of date, i.e. there is a newer version available through conda.")
text.append(f"I have updated {args.repo} {update}.")

with open(args.shed) as f:
    y = yaml.load(f, Loader=yaml.SafeLoader)

if y.get('homepage_url'):
    url = y.get('homepage_url').strip('/')
    if 'github.com' in url:
        if len(url.split('github.com')[1].split('/')) > 1:
            url += '/releases'
    text.append(f'**Project home page:** {url}')

if y.get('maintainers'):
    text.append('**Maintainers:** ' + ', '.join([f'@{m}' for m in y.get('maintainers')]))

text.append("For any comments, queries or criticism about the bot, not related to the tool being updated in this PR, please create an issue [here](https://github.com/galaxyproject/tools-iuc/new).")

with open(args.out, 'w') as f:
    f.write('\n\n'.join(text))

print(f'Updating {args.repo} {update}')
