import pytest

from tssenrich import (
    generate_tss, generate_flanks, flanks_bed_str, flanks_bed_tool
)

@pytest.fixture()
def flanks():
    return flanks_bed_tool(flanks_bed_str(generate_flanks(generate_tss())))


def test_flanks(flanks):
    flanks.head()
