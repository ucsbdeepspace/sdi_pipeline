"""initial test

Revision ID: 96d42113364c
Revises: 
Create Date: 2020-09-08 17:38:53.105003

"""
from alembic import op
import sqlalchemy as sa


# revision identifiers, used by Alembic.
revision = '96d42113364c'
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    op.add_column('Image', sa.Column('coeff_a', sa.Float))
    op.add_column('Image', sa.Column('coeff_b', sa.Float))


def downgrade():
    op.drop_column('Image', 'coeff_a')
    op.drop_column('Image', 'coeff_b')
