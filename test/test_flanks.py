import pytest

from tssenrich import (
    generate_tss, generate_tss_flanks, tss_flanks_bed_str, tss_flanks_bed_tool
)

@pytest.fixture()
def flanks():
    return tss_flanks_bed_tool(
        tss_flanks_bed_str(generate_tss_flanks(generate_tss()))
    )


def test_flanks(flanks):
    flanks.head()
